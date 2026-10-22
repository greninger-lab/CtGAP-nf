/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTP - Read trimming and quality control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: fastp (1-trim.smk)
    
    Performs adapter detection (PE mode), sliding-window quality trimming
    from both ends, and right-side cutting.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTP {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/fastp:1.4.0--h43da1c4_0"

    publishDir "${params.outdir}/${sample_id}/trim", mode: 'copy', pattern: '*_trim_R{1,2}.fastq.gz'
    publishDir "${params.outdir}/${sample_id}/trim", mode: 'copy', pattern: '*.{html,json}'
    publishDir "${params.outdir}/status",            mode: 'copy', pattern: '*.status.txt'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_trim_R1.fastq.gz"), path("${sample_id}_trim_R2.fastq.gz"), emit: reads
    path "${sample_id}.html",                                                                            emit: html
    path "${sample_id}.json",                                                                            emit: json
    path "trim.${sample_id}.status.txt",                                                                 emit: status

    script:
    """
    fastp \\
        -i ${r1} \\
        -I ${r2} \\
        -o ${sample_id}_trim_R1.fastq.gz \\
        -O ${sample_id}_trim_R2.fastq.gz \\
        -h ${sample_id}.html \\
        -j ${sample_id}.json \\
        --detect_adapter_for_pe \\
        --cut_front \\
        --cut_front_window_size=${params.fastp_window_size} \\
        --cut_front_mean_quality=${params.fastp_mean_quality} \\
        --cut_tail \\
        --cut_tail_window_size=${params.fastp_window_size} \\
        --cut_tail_mean_quality=${params.fastp_mean_quality} \\
        --cut_right \\
        --cut_right_window_size=${params.fastp_right_window_size} \\
        --cut_right_mean_quality=${params.fastp_right_mean_quality} \\
        --thread ${task.cpus}

    touch trim.${sample_id}.status.txt
    """
}
