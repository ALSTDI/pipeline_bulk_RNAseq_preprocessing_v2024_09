process checksum_compare {
    /* process should not give an error,
    even when the files are not identical,
    because cmp command is inside the if statement */
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    cpus 4
    memory '16 GB'
    publishDir params.out_bucket, mode: 'copy'
    
    input:
    tuple val(SAMPLE), path(SAMPLE_DIR), path(SAMPLE_MD5)

    output:
    tuple val(SAMPLE), 
    env(is_faithful), emit: is_faithful

    shell:
    """
    md5sum !{SAMPLE_DIR}/*.gz > test_sum
    cat test_sum | tr ' ' '\t' | cut -f 1 | sort -g > file1
    cat !{SAMPLE_MD5} |  tr ' ' '\t' | cut -f 1 | sort -g > file2
    
    if cmp --silent file1 file2;
    then
        echo !{SAMPLE} is good to go > !{SAMPLE}_md5.summary
        is_faithful=true
    else
        echo !{SAMPLE} files are CORRUPTED > !{SAMPLE}_md5.summary
        is_faithful=false
        
    fi 
    """
}

process fq_line_compare_paired {
    /* process should not give an error,
    even when the files are not identical,
    because cmp command is inside the if statement */
    /* use this instead of checksum_compare in the case when
    there is no available historic checksums to compare to */
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    cpus 4
    memory '16 GB'
    
    input:
    tuple val(SAMPLE), path(R1), path(R2)

    output:
    tuple val(SAMPLE), 
    env(is_faithful), emit: is_faithful

    shell:
    """
    zcat !{R1} > !{SAMPLE}_R1.fastq
    zcat !{R2} > !{SAMPLE}_R2.fastq

    wc -l !{SAMPLE}_R1.fastq | awk '{print \$1}' > !{SAMPLE}_r1_reads
    wc -l !{SAMPLE}_R2.fastq | awk '{print \$1}' > !{SAMPLE}_r2_reads

    read_token=\$(wc -l !{SAMPLE}_R1.fastq | awk '{print \$1}')

    if cmp --silent !{SAMPLE}_r1_reads !{SAMPLE}_r2_reads;
    then
            if [ "\$(( \$read_token % 4 ))" -eq 0 ]; then
                echo !{SAMPLE} has the same number of reads between R1 and R2 > !{SAMPLE}_md5.summary
                echo !{SAMPLE} reads are divisible by 4 >> !{SAMPLE}_md5.summary
                is_faithful=true
            else
                echo !{SAMPLE} files are CORRUPTED > !{SAMPLE}_md5.summary
                is_faithful=false
            fi
    else
        echo !{SAMPLE} files are CORRUPTED > !{SAMPLE}_md5.summary
        is_faithful=false
    fi 
    """
}

process fq_line_compare_single {
    /* use this instead of checksum_compare in the case when
    there is no available historic checksums to compare to 
    for SINGLE-END reads */
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'
    cpus 4
    memory '16 GB'
    
    input:
    tuple val(SAMPLE), path(R1)

    output:
    tuple val(SAMPLE), 
    env(is_faithful), emit: is_faithful

    shell:
    """
    zcat !{R1} > !{SAMPLE}_R1.fastq

    read_token=\$(wc -l !{SAMPLE}_R1.fastq | awk '{print \$1}')
    if [ "\$(( \$read_token % 4 ))" -eq 0 ]; then
        echo !{SAMPLE} reads are divisible by 4 > !{SAMPLE}_md5.summary
        is_faithful=true
    else
        echo !{SAMPLE} files are CORRUPTED > !{SAMPLE}_md5.summary
        is_faithful=false
    fi
    """
}

process sample_fasta_paired {
    /* Subset by default 100k reads for faster run with salmon to detect library types */
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/seqtk_compress'
    cpus 4
    memory '16 GB'

    input: 
    tuple val(SAMPLE), path(R1), path(R2)
    path trim_R1
    path trim_R2
    val NR_READS

    output: 
    tuple val("${SAMPLE}"), path("${SAMPLE}_sub_R1.fastq.gz"), path("${SAMPLE}_sub_R2.fastq.gz")
    script:
    """
    seqtk sample -s2358 ${trim_R1} ${NR_READS} | gzip > ${SAMPLE}_sub_R1.fastq.gz
    seqtk sample -s2358 ${trim_R2} ${NR_READS} | gzip > ${SAMPLE}_sub_R2.fastq.gz
    """
}

process sample_fasta_single {
    /* Subset by default 100k reads for faster run with salmon to detect library types */
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/seqtk_compress'
    cpus 4
    memory '16 GB'

    input: 
    tuple val(SAMPLE), path(R1)
    path trim_R1
    val NR_READS

    output: 
    tuple val("${SAMPLE}"), path("${SAMPLE}_sub_R1.fastq.gz")

    script:
    """
    seqtk sample -s2358 ${trim_R1} ${NR_READS} | gzip > ${SAMPLE}_sub_R1.fastq.gz
    """
}

process MULTIQC {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/multiqc'
    cpus 4
    memory '16 GB'
    publishDir "${params.out_bucket}", mode: 'copy'

    input:
    path '*'
    val pasi

    output:
    path "multiqc_report_${pasi}.html"

    script:
    """
    multiqc .
    mv multiqc_report.html multiqc_report_${pasi}.html
    """
}