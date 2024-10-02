include { checksum_compare; sample_fasta_single; sample_fasta_paired } from './modules/diag'
include { fastp_paired; fastqc_paired; STAR_paired; salmon_paired; rsem_expr_paired }  from './modules/paired_processes'
include { fastp_single; fastqc_single; STAR_single; salmon_single; rsem_expr_single }  from './modules/single_processes'
include { MULTIQC as MULTIP; MULTIQC as MULTIS } from './modules/diag'
/*
Workflow seeks to do the following:
1. Determine if the fastq files are faithful (i.e. not corrupted) with checksum_compare
2. Failed samples will be written out to sample_failed.csv with collectFile()
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
workflow CHECK_MD5 {
    main:
    input_ch=INPUT().map{
        sample, sampleDir, sampleChecksum, R1, R2 ->
        [ sample, sampleDir, sampleChecksum ]
    }

    // TODO: https://github.com/nextflow-io/nextflow/issues/4446
    // Implemented after 24.05 but we don't know how to use it yet. 
    checksum_compare(input_ch)

    emit:
    check_res=checksum_compare.out
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

workflow {
    input_ch = INPUT().map{ sample_name, sample_dir, sample_checksum, read1, read2 
                            ->
                            [ sample_name, read1, read2] }
    
    pasi_ch = PAIRED_OR_SINGLE()

    reads_ch = pasi_ch.join(input_ch, remainder: true)
    
    CHECK_MD5()
    // Sort checksum outputs into faithful samples and vice versa - with the branch operator
    CHECK_MD5.out.check_res \
    .branch {
        failed: it[1] == "false"
        passed: it[1] == "true"
    }.set{ result }
    
    // Write corrupted samples out for traceback later
    result.failed.map { sample, is_faithful -> sample }
                 .collectFile(name: "sample_corrupted.txt", storeDir: '.' , newLine: true )
    
    // Use the join operator to make a new input channel, with only the samples that passed the checksum test 
    // and their respective reads 
    // input_ch.view()
    // reads_ch.view()
    meta_ch = result.passed.join(reads_ch, remainder: false)

    meta_ch.branch{
        paired: it[2] == "paired"
        single: it[2] == "single"
    }.set{ result2 }

    paired_ch = result2.paired.map{ sample_name, pass, paired, read1, read2 -> [sample_name, read1, read2] }
    single_ch = result2.single.map{ sample_name, pass, paired, read1, read2 -> [sample_name, read1] }

    /* FINALLY EVOKE the PAIRED-END WORKFLOW */
    // # Trim reads
    fastp_paired(paired_ch)
    // # Subset reads and predict strandedness
    sample_fasta_paired( paired_ch, fastp_paired.out.trim_r1, fastp_paired.out.trim_r2, 50000 )
    salmon_paired( Channel.fromPath(params.salmonDir), sample_fasta_paired.out )
    // # Alignment and expression quantification
    STAR_paired( paired_ch, Channel.fromPath(params.refDir), fastp_paired.out.trim_r1, fastp_paired.out.trim_r2 )
    rsem_expr_paired( Channel.fromPath(params.refDir), STAR_paired.out.Tr_bam, salmon_paired.out.salmon_out  )
    // QC report
    fastqc_paired( paired_ch, fastp_paired.out.trim_r1, fastp_paired.out.trim_r2 )

    /* SIMILAR SINGLE-END WORKFLOW */
    fastp_single(single_ch)
    sample_fasta_single(single_ch, fastp_single.out.trim_r1, 50000)
    salmon_single( Channel.fromPath(params.salmonDir), sample_fasta_single.out )
    STAR_single( single_ch, Channel.fromPath(params.refDir), fastp_single.out.trim_r1 )
    rsem_expr_single( Channel.fromPath(params.refDir), STAR_single.out.Tr_bam, salmon_single.out.salmon_out )
    fastqc_single_ch = fastqc_single(single_ch, fastp_single.out.trim_r1)

    //TODO: multiQC
    MULTIP(fastqc_paired.out[0]
            .mix(fastqc_paired.out[1])
            .mix(STAR_paired.out[0])
            .mix(rsem_expr_paired.out)
            .collect(), "paired")

    MULTIS(fastqc_single.out[0]
            .mix(fastqc_single.out[1])
            .mix(STAR_single.out[0])
            .mix(rsem_expr_single.out)
            .collect(), "single")
}