/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_BOWTIE - Align scrubbed reads to specified reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: bowtie (5-align.smk)
    Split into two processes for container isolation:
      1. REF_BOWTIE_ALIGN    - bowtie2 alignment (outputs SAM)
      2. REF_SAMTOOLS_SORT    - samtools name-sort + flagstat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_BOWTIE_ALIGN {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.ref.sam"), emit: sam

    script:
    """
    bowtie2 -x ${params.ref_genome_bt2_prefix} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --threads ${task.cpus} \\
        -S ${sample_id}.ref.sam \\
        2> bowtie2.log
    """
}

process REF_SAMTOOLS_SORT {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/bowtie2", mode: 'copy'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam"),          emit: bam
    tuple val(sample_id), path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    """
    samtools view -bSh ${sam} | \\
        samtools sort -n -@${task.cpus} \\
        -o ${sample_id}.bam

    samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat.txt

    # Check if any reads mapped
    mapped=\$(grep "mapped (" ${sample_id}.flagstat.txt | head -1 | awk '{print \$1}')
    if [ "\${mapped:-0}" -eq 0 ]; then
        echo "WARNING: No reads mapped to reference for ${sample_id}" >> ${sample_id}.flagstat.txt
    fi
    """
}
