/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KRAKEN2_CLASSIFY - Taxonomic classification of host-depleted reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: kraken2_classify (2-scrub_hocort.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process KRAKEN2_CLASSIFY {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/kraken2:2.17.2--pl5321h3be2455_0"

    publishDir "${params.outdir}/${sample_id}/scrub", mode: 'copy', pattern: '*.kraken2.report'
    publishDir "${params.outdir}/${sample_id}/scrub", mode: 'copy', pattern: '*.kraken2.labels'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.kraken2"), path("${sample_id}.kraken2.report"), emit: results
    path "${sample_id}.kraken2.labels",                                                       emit: labels

    script:
    """
    kraken2 \\
        --db ${params.kraken_db} \\
        --threads ${task.cpus} \\
        --report ${sample_id}.kraken2.report \\
        --paired ${r1} ${r2} \\
        --output ${sample_id}.kraken2 \\
        --use-names

    # Compact labels for KrakenTools (handle empty output)
    if [ -s ${sample_id}.kraken2 ]; then
        cut -f2,3 ${sample_id}.kraken2 > ${sample_id}.kraken2.labels
    else
        touch ${sample_id}.kraken2.labels
    fi
    """
}
