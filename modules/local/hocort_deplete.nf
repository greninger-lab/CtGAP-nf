/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HOCORT_DEPLETE - Host read depletion via HoCoRT (bowtie2 backend)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: host_deplete_hocort (2-scrub_hocort.smk)
    Keeps UNMAPPED reads (non-host).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HOCORT_DEPLETE {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/hocort:1.2.2--py39hdfd78af_0"

    publishDir "${params.outdir}/${sample_id}/scrub", mode: 'copy', pattern: '*_hostdepleted_R{1,2}.fastq.gz'
    publishDir "${params.outdir}/status",             mode: 'copy', pattern: '*.status.txt'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_hostdepleted_R1.fastq.gz"), path("${sample_id}_hostdepleted_R2.fastq.gz"), emit: reads
    path "hocort.${sample_id}.status.txt",                                                                               emit: status

    script:
    def extra_args = params.hocort_bowtie2_args ? "-c=\"${params.hocort_bowtie2_args}\"" : ""
    def opts       = params.hocort_opts ?: ""
    """
    hocort map bowtie2 \\
        -x ${params.host_bt2_prefix} \\
        -i ${r1} ${r2} \\
        -o ${sample_id}_hostdepleted_R1.fastq.gz ${sample_id}_hostdepleted_R2.fastq.gz \\
        -t ${task.cpus} \\
        -p end-to-end \\
        -f true \\
        ${extra_args} \\
        ${opts}

    touch hocort.${sample_id}.status.txt
    """
}
