/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTP_QC - Post-scrub quality control report (no trimming)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: fastp_qc (3-qc.smk)
    Runs fastp in report-only mode on scrubbed reads to generate QC metrics.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTP_QC {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/fastp:1.4.0--h43da1c4_0"

    publishDir "${params.outdir}/${sample_id}/qc",     mode: 'copy', pattern: '*.{html,json}'
    publishDir "${params.outdir}/${sample_id}/status",  mode: 'copy', pattern: '*.status.txt'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_scrub.html"), path("${sample_id}_scrub.json"), emit: reports
    path "${sample_id}_scrub.json",                                                          emit: json
    path "fastp_qc.${sample_id}.status.txt",                                                 emit: status

    script:
    """
    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -h ${sample_id}_scrub.html \\
        -j ${sample_id}_scrub.json \\
        --disable_adapter_trimming \\
        --disable_quality_filtering \\
        --disable_length_filtering \\
        --thread ${task.cpus}

    touch fastp_qc.${sample_id}.status.txt
    """
}
