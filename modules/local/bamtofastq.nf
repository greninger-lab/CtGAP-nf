/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAMTOFASTQ - Extract aligned reads from reference-guided alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: bamtofastq (5-align.smk)
    Converts name-sorted BAM back to paired FASTQ for re-assembly.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BAMTOFASTQ {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/bedtools:2.27.1--h077b44d_9"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/bamtofastq", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_refg_r1.fastq.gz"), path("${sample_id}_refg_r2.fastq.gz"), emit: reads

    script:
    """
    bedtools bamtofastq \\
        -i ${bam} \\
        -fq ${sample_id}_refg_r1.fastq \\
        -fq2 ${sample_id}_refg_r2.fastq

    gzip -f ${sample_id}_refg_r1.fastq
    gzip -f ${sample_id}_refg_r2.fastq
    """
}
