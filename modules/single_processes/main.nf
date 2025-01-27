process fastp_single {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    
    publishDir "${params.out_bucket}/${SAMPLE}/logs/", mode: 'copy', pattern: '*.json'

    input:
    tuple val(SAMPLE), path(R1)

    output:
    path("${SAMPLE}_fastp.json"), emit: fastp_log
    path("clt_${SAMPLE}_R1.fastq.gz"), emit: trim_r1

    script:
    """
    fastp \
            --in1 ${R1} \
            --out1 clt_${SAMPLE}_R1.fastq.gz \
            --length_required 35 \
            --cut_front_window_size 1 \
            --cut_front_mean_quality 13 \
            --cut_front \
            --cut_tail_window_size 1 \
            --cut_tail_mean_quality 13 \
            --cut_tail \
            -w 2 -y -r \
            -j ${SAMPLE}_fastp.json
    """
}



process fastqc_single {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    
    publishDir "${params.out_bucket}/${SAMPLE}/logs/", mode: 'copy'

    input:
    tuple val(SAMPLE), path(R1)
    path(trim_R1)

    output:
    path("raw_fastqc_${SAMPLE}")
    path("clean_fastqc_${SAMPLE}")

    script:
    """
    mkdir raw_fastqc_${SAMPLE}
    fastqc \
                        -o raw_fastqc_${SAMPLE} \
                        -f fastq \
                        -q \
                        -t ${task.cpus} \
                        ${R1} ;

    mkdir clean_fastqc_${SAMPLE}
    fastqc \
                        -o clean_fastqc_${SAMPLE} \
                        -f fastq \
                        -q \
                        -t ${task.cpus} \
                        ${trim_R1} ;
    """
}

process STAR_single {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    
    publishDir "${params.out_bucket}/${SAMPLE}/logs/STAR/", mode: 'copy', pattern: '*.{out,log}'

    input: 
    tuple val(SAMPLE), path(R1)
    path STARREF
    path trim_R1

    output: 
    path("*")
    path("ctl_${SAMPLE}Aligned.toTranscriptome.out.bam"), emit: Tr_bam

    script: 
    """
    STAR \
                        --twopassMode Basic \
                        --readFilesCommand zcat \
                        --runThreadN ${task.cpus}-2 \
                        --runMode alignReads \
                        --genomeDir ${STARREF} \
                        --readFilesIn ${trim_R1} \
                        --outSAMtype BAM Unsorted \
                        --outSAMstrandField intronMotif \
                        --outFileNamePrefix ctl_${SAMPLE} \
                        --outFilterIntronMotifs RemoveNoncanonical \
                        --outReadsUnmapped Fastx \
			            --quantMode TranscriptomeSAM ;
    """
}

process salmon_single {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    
    publishDir "${params.out_bucket}/${SAMPLE}/logs", mode: 'copy'

    input:
    path salmon_index
    tuple val(SAMPLE), path(sub_R1)

    output: 
    path("salmon_quant.log")
    tuple val(SAMPLE), env(rsemSTRAND), emit: salmon_out
    

    shell: 
    '''
    salmon quant \
                        -i !{salmon_index} \
                        -l A \
                        --seqBias \
                        --gcBias \
                        --posBias \
                        -r !{sub_R1} \
                        -p !{task.cpus} \
                        -o salmon_!{SAMPLE}

    mv salmon_!{SAMPLE}/logs/salmon_quant.log salmon_quant.log
    # Define strandedness by salmon output
    salmonOUTSTRAND="$(grep Automatically salmon_quant.log | tr ' ' '\t' | cut -f 12)"

    if [ "${salmonOUTSTRAND}" == "SF" ]
    then
        rsemSTRAND="forward"
    elif [ "${salmonOUTSTRAND}" == "SR" ]
    then
        rsemSTRAND="reverse"
    elif [ "${salmonOUTSTRAND}" == "U" ]
    then
        rsemSTRAND="none"
    fi

    '''
}

process rsem_expr_single {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'

    publishDir(path: {"${params.out_bucket}/${SAMPLE}/RSEM"}, mode: 'copy', pattern: '*.results')
    publishDir(path: {"${params.out_bucket}/${SAMPLE}/logs"}, mode: 'copy', pattern: '*.log')
    publishDir(path: {"${params.out_bucket}/${SAMPLE}/RSEM"}, mode: 'copy', pattern: '*.out.bam')

    input:
    path(rsem_path)
    path(Tr_bam)
    tuple val(SAMPLE), val(rsemSTRAND)

    output: 
    path("*")

    shell: 
    '''
    # Prepare STAR bam file 
    samtools view -H !{Tr_bam} > file1 ;
    samtools view -@ !{task.cpus} !{Tr_bam} | awk '{printf "%s", \$0 " "; getline; print}' | sort -S 20G -T ./ | tr ' ' '\n' > file2 ;
    cat file1 file2 | samtools view -@ !{task.cpus} -bS - > rsem_!{SAMPLE}_Aligned.toTranscriptome.out.bam ;

    # RSEM index is supplied weirdly as path_to_index/index_prefix (*.grp, *.chrlist, *idx.fa, etc.)
    # Learned from nf-core rna-seq how to supply rsem index
    INDEX="$(find -L ./ -name "*.grp" | sed 's/\\.grp\$//')"
    echo ${INDEX}

    # Finally initiate rsem
    rsem-calculate-expression \
                        --strandedness !{rsemSTRAND} \
                        --no-bam-output \
                        --estimate-rspd \
                        --alignments \
                        rsem_!{SAMPLE}_Aligned.toTranscriptome.out.bam \
                        ${INDEX} \
                        !{SAMPLE} >& rsem_!{SAMPLE}.log

    sed -i '1i Quantifying expression level for !{SAMPLE} with strandedness set to !{rsemSTRAND}' rsem_!{SAMPLE}.log
    '''
}
