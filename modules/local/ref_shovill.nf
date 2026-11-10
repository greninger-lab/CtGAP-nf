/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_SHOVILL - De novo assembly of reference-guided extracted reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_shovill (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_SHOVILL {
    tag "${sample_id}"
    label 'process_high'

    container "quay.io/biocontainers/shovill:1.4.2--hdfd78af_1"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/shovill", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("contigs.fa"), emit: contigs
    path "shovill.log",                       emit: log

    script:
    """
    # Preserve Shovill's exit status when logging through tee.
    set -o pipefail

    shovill \\
        --R1 ${r1} \\
        --R2 ${r2} \\
        --outdir shovill_tmp \\
        --gsize ${params.shovill_gsize} \\
        --depth ${params.shovill_depth} \\
        --force \\
        --cpus ${task.cpus} 2>&1 | tee shovill.log

    if [ ! -s shovill_tmp/contigs.fa ]; then
        printf '%s\\n' 'REF_SHOVILL: no non-empty contigs.fa was produced for ${sample_id}; inspect shovill.log and the input read counts.' >&2
        exit 1
    fi

    mv shovill_tmp/contigs.fa .
    """
}
