include { sample_fasta_single; sample_fasta_paired } from './modules/diag'
include { fq_line_compare_paired; fq_line_compare_single } from './modules/diag'
include { fastp_paired; fastqc_paired; STAR_paired; salmon_paired; rsem_expr_paired }  from './modules/paired_processes'
include { fastp_single; fastqc_single; STAR_single; salmon_single; rsem_expr_single }  from './modules/single_processes'
include { MULTIQC as MULTIP; MULTIQC as MULTIS } from './modules/diag'
/*
Workflow seeks to do the following:
1. Determine if the fastq files are faithful (i.e. not corrupted) with fq_line_compare
2. Failed samples will be written out to "sample_corrupted_*.txt" with collectFile()
3. Passed samples will be passed on to evaluate if a sample is paired- or single-end
4. Emit different process workflow based on this information. 
*/
workflow INPUT {
    main: 
    input_ch=Channel.fromPath("sampleSheet.csv").splitCsv(header: true)
    | map {
        row -> [ row.sample_name, //tuple
                 file("${params.projectDir}/${row.sample_dir}"),
                 file("${params.projectDir}/${row.sample_dir}/${row.checksum}"),
                 file("${params.projectDir}/${row.sample_dir}/${row.r1}"),
                 file("${params.projectDir}/${row.sample_dir}/${row.r2}") ]
          }
    emit:
    input_ch
}

workflow PAIRED_OR_SINGLE {
    /* Evaluate whether this is a paired- or single-end experiment 
    Based on USER-SUPPLED sampleSheet.csv 
    If there is not R2, it's single-end
    The reason why we do not want a "paired_or_single" process is because sometimes file naming convention is
    not _R2 but _2 - and this will require hard-corded elements to be changed in the workflow 
    As such, we leave it to the users to determine this 
    */
    main:
    // Define R1 and R2 outputs
    reads_ch = Channel.fromPath("sampleSheet.csv").splitCsv( header: true ) 
    | map {
        row ->
        // This is the csv entry
        // So we don't need the full path
        // Just want to use this to determine whether paired or single 
        r1=row.r1
        r2=row.r2
        // conditional statement 
        is_paired =  r2.toString()=='' ? 'single' : 'paired'
        // Return sample and paired or single 
        [row.'sample_name', is_paired]
        }
    
    emit:
    reads_ch
    
}

workflow CHECK_INTEGRITY_BY_LINES {
    main:
    input_ch=INPUT().map{
        sample, sampleDir, sampleChecksum, R1, R2 ->
        [ sample, sampleDir, sampleChecksum, R1, R2 ]
    }

    pasi_ch = PAIRED_OR_SINGLE()

    reads_ch = pasi_ch.join(input_ch, remainder: true)

    reads_ch.branch{
        paired: it[1] == "paired"
        single: it[1] == "single"
    }.set{ result2 }

    paired_ch = result2.paired.map{ sample_name, paired, sample_dir, checksum, read1, read2 -> [sample_name, read1, read2] }
    single_ch = result2.single.map{ sample_name, paired, sample_dir, checksum, read1, read2 -> [sample_name, read1] }

    // CHECK INTEGRITY OF FASTQS by the NUMBER of LINES
    fq_line_compare_paired(paired_ch)
    fq_line_compare_single(single_ch)

    emit:
    check_res_paired=fq_line_compare_paired.out
    check_res_single=fq_line_compare_single.out   
}

workflow PAIRED {
    input_ch = INPUT().map{ sample_name, sample_dir, sample_checksum, read1, read2 
                            ->
                            [ sample_name, read1, read2] }
    
    CHECK_INTEGRITY_BY_LINES()
    // Sort checksum outputs into faithful samples and vice versa - with the branch operator
    CHECK_INTEGRITY_BY_LINES.out.check_res_paired \
    .branch {
        failed: it[1] == "false"
        passed: it[1] == "true"
    }.set{ result }
    
    // Write corrupted samples out for traceback later
    result.failed.map { sample, is_faithful -> sample }
                 .collectFile(name: "sample_corrupted_paired.txt", storeDir: '.' , newLine: true )
    
    // Use the join operator to make a new input channel, with only the paired-end samples that passed the checksum test 
    // and their respective reads 
    paired_ch = result.passed.join(input_ch, remainder: false).map{
        sample_name, pass, read1, read2 -> [sample_name, read1, read2]
    }
    
    /* FINALLY EVOKE the PAIRED-END WORKFLOW */
    // # Trim reads
    fastp_paired(paired_ch)
    // # Subset reads and predict strandedness
    sample_fasta_paired( paired_ch, fastp_paired.out.trim_r1, fastp_paired.out.trim_r2, 50000 )
    salmon_paired( file(params.salmonDir), sample_fasta_paired.out )
    // # Alignment and expression quantification
    STAR_paired( paired_ch, file(params.refDir), fastp_paired.out.trim_r1, fastp_paired.out.trim_r2 )
    rsem_expr_paired( file(params.refDir), STAR_paired.out.Tr_bam, salmon_paired.out.salmon_out  )
    // QC report
    fastqc_paired( paired_ch, fastp_paired.out.trim_r1, fastp_paired.out.trim_r2 )

    //MultiQC Report
    MULTIP(fastqc_paired.out[0]
            .mix(fastqc_paired.out[1])
            .mix(STAR_paired.out[0])
            .mix(rsem_expr_paired.out)
            .collect(), "paired")
}

workflow SINGLE {
                    /* SIMILAR SINGLE-END WORKFLOW */
    input_ch = INPUT().map{ sample_name, sample_dir, sample_checksum, read1, read2 
                            ->
                            [ sample_name, read1] }
    
    CHECK_INTEGRITY_BY_LINES()
    // Sort checksum outputs into faithful samples and vice versa - with the branch operator
    CHECK_INTEGRITY_BY_LINES.out.check_res_single \
    .branch {
        failed: it[1] == "false"
        passed: it[1] == "true"
    }.set{ result }
    
    // Write corrupted samples out for traceback later
    result.failed.map { sample, is_faithful -> sample }
                 .collectFile(name: "sample_corrupted_single.txt", storeDir: '.' , newLine: true )
    
    // Use the join operator to make a new input channel, with only the SINGLE-end samples that passed the checksum test 
    // and their respective reads 
    single_ch = result.passed.join(input_ch, remainder: false).map{
        sample_name, pass, read1 -> [sample_name, read1]
    }
    
                /* EVOKE QUANTIFICATION WORKFLOW */
    fastp_single(single_ch)
    sample_fasta_single(single_ch, fastp_single.out.trim_r1, 50000)
    salmon_single( file(params.salmonDir), sample_fasta_single.out )
    STAR_single( single_ch, file(params.refDir), fastp_single.out.trim_r1 )
    rsem_expr_single( file(params.refDir), STAR_single.out.Tr_bam, salmon_single.out.salmon_out )

    fastqc_single_ch = fastqc_single(single_ch, fastp_single.out.trim_r1)
    MULTIS(fastqc_single.out[0]
            .mix(fastqc_single.out[1])
            .mix(STAR_single.out[0])
            .mix(rsem_expr_single.out)
            .collect(), "single")
}

workflow {
    PAIRED()
    SINGLE()
}