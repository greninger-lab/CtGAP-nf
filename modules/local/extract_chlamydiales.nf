/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT_CHLAMYDIALES - Extract Chlamydiales reads using KrakenTools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: extract_chlamydiales_reads (2-scrub_hocort.smk)
    Extracts reads classified as Chlamydiales (taxid 51291) + children + parents.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EXTRACT_CHLAMYDIALES {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/jefffurlong/krakentools"

    publishDir "${params.outdir}/${sample_id}/scrub", mode: 'copy', pattern: '*_trim_scrub_R{1,2}.fastq.gz'
    publishDir "${params.outdir}/status",             mode: 'copy', pattern: '*.status.txt'

    input:
    tuple val(sample_id), path(r1), path(r2), path(k2_output), path(k2_report)

    output:
    tuple val(sample_id), path("${sample_id}_trim_scrub_R1.fastq.gz"), path("${sample_id}_trim_scrub_R2.fastq.gz"), emit: reads
    path "scrub.${sample_id}.status.txt",                                                                            emit: status

    script:
    """
    tmp1=\$(mktemp)
    tmp2=\$(mktemp)

    python3 /usr/share/KrakenTools/extract_kraken_reads.py \\
        -k ${k2_output} \\
        --report ${k2_report} \\
        -s ${r1} -s2 ${r2} \\
        -o "\$tmp1" -o2 "\$tmp2" \\
        -t ${params.chlamydiales_taxid} \\
        --include-children \\
        --include-parents \\
        --fastq-output

    gzip -c "\$tmp1" > ${sample_id}_trim_scrub_R1.fastq.gz
    gzip -c "\$tmp2" > ${sample_id}_trim_scrub_R2.fastq.gz
    rm -f "\$tmp1" "\$tmp2"

    touch scrub.${sample_id}.status.txt
    """
}
