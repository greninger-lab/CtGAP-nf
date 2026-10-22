/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC - Aggregate QC reports across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: multiqc (3-qc.smk)
    Collects all fastp JSON reports and produces a combined MultiQC report.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MULTIQC {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path ('fastp_reports/*')

    output:
    path "multiqc_report.html",      emit: report
    path "multiqc_report_data",      emit: data

    script:
    """
    multiqc fastp_reports/ \\
        -o . \\
        -n multiqc_report
    """
}
