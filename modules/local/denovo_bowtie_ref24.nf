/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DENOVO_BOWTIE_REF24 - Align scrubbed reads against reference set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: denovo_bowtie_ref24 (4-assemble.smk)
    Split into two processes for container isolation:
      1. DENOVO_BOWTIE_ALIGN  - bowtie2 alignment (outputs SAM)
      2. DENOVO_SAMTOOLS_REF24 - samtools sort/index/coverage/flagstat + awk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DENOVO_BOWTIE_ALIGN {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.ref24.sam"), emit: sam

    script:
    """
    bowtie2 -x ${params.refset_bt2_prefix} \\
        -1 ${r1} \\
        -2 ${r2} \\
        -p ${task.cpus} \\
        -S ${sample_id}.ref24.sam \\
        2> bowtie2.log
    """
}

process DENOVO_SAMTOOLS_REF24 {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

    publishDir "${params.outdir}/${sample_id}/denovo/bowtie_ref24", mode: 'copy'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.ref24.bam"), path("${sample_id}.ref24.bam.bai"),  emit: bam
    tuple val(sample_id), path("${sample_id}.coverage.ref24.tsv"),                              emit: coverage
    tuple val(sample_id), path("${sample_id}.flagstat.txt"),                                    emit: flagstat

    script:
    """
    samtools view -bSh ${sam} | \\
        samtools sort -@${task.cpus} -o ${sample_id}.ref24.bam

    samtools index ${sample_id}.ref24.bam

    samtools coverage ${sample_id}.ref24.bam > tmp.cov.tsv

    # Add sample_name column; handle empty coverage (header-only) gracefully
    if [ \$(wc -l < tmp.cov.tsv) -gt 1 ]; then
        awk -v sample="${sample_id}" 'BEGIN{OFS="\\t"} NR==1{print "sample_name", \$0} NR>1{print sample, \$0}' tmp.cov.tsv > ${sample_id}.coverage.ref24.tsv
    else
        echo -e "sample_name\\t\$(head -1 tmp.cov.tsv)" > ${sample_id}.coverage.ref24.tsv
    fi

    samtools flagstat ${sample_id}.ref24.bam > ${sample_id}.flagstat.txt
    """
}
