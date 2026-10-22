/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GAPFILLER - Fill gaps in scaffolded assembly using Gap2Seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: gapfiller (4-assemble.smk)
    Falls back to scaffolded assembly if Gap2Seq fails.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GAPFILLER {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/gap2seq:3.1.1a--py311hc84137b_6"

    publishDir "${params.outdir}/${sample_id}/denovo/gap2seq", mode: 'copy', pattern: 'filled.fasta'

    input:
    tuple val(sample_id), path(scaffold), path(r1), path(r2)

    output:
    tuple val(sample_id), path("filled.fasta"), emit: filled

    script:
    """
    set +e
    Gap2Seq \\
        --scaffolds ${scaffold} \\
        --filled filled.fasta \\
        --reads ${r1},${r2} \\
        --threads ${task.cpus} > gap2seq.log 2>&1
    exitcode=\$?
    set -e

    if [ \$exitcode -ne 0 ]; then
        cp ${scaffold} filled.fasta
    fi
    """
}
