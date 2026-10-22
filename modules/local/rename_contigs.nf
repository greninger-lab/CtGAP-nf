/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RENAME_CONTIGS - Rename contigs with sample-prefixed names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: rename (4-assemble.smk)
    Uses seqkit to rename contigs as {sample}_contig1, {sample}_contig2, etc.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RENAME_CONTIGS {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/seqkit:2.14.0--hb192632_0"

    publishDir "${params.outdir}/${sample_id}/denovo", mode: 'copy'

    input:
    tuple val(sample_id), path(filled)

    output:
    tuple val(sample_id), path("denovo_${sample_id}.final.fasta"), emit: assembly

    script:
    """
    seqkit replace -p "(.*)" -r '${sample_id}_contig{nr}' ${filled} > denovo_${sample_id}.final.fasta
    """
}
