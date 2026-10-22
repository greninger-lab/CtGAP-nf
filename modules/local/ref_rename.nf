/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_RENAME - Rename contigs for reference-guided assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_rename (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_RENAME {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/seqkit:2.14.0--hb192632_0"

    publishDir "${params.outdir}/${sample_id}/ref-denovo", mode: 'copy'

    input:
    tuple val(sample_id), path(filled)

    output:
    tuple val(sample_id), path("ref-denovo_${sample_id}.final.fasta"), emit: assembly

    script:
    """
    seqkit replace -p "(.*)" -r '${sample_id}_contig{nr}' ${filled} > ref-denovo_${sample_id}.final.fasta
    """
}
