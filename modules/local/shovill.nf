/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SHOVILL - De novo assembly using Shovill (SPAdes-based)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: shovill (4-assemble.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SHOVILL {
    tag "${sample_id}"
    label 'process_high'

    container "quay.io/biocontainers/shovill:1.4.2--hdfd78af_1"

    publishDir "${params.outdir}/${sample_id}/denovo/shovill", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("contigs.fa"),   emit: contigs
    path "shovill.log",                         emit: log

    script:
    """
    shovill \\
        --R1 ${r1} \\
        --R2 ${r2} \\
        --outdir shovill_tmp \\
        --gsize ${params.shovill_gsize} \\
        --depth ${params.shovill_depth} \\
        --force \\
        --cpus ${task.cpus} 2>&1 | tee shovill.log

    # Move outputs to working dir for publishDir
    mv shovill_tmp/contigs.fa .
    """
}
