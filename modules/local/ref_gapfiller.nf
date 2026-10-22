/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_GAPFILLER - Fill gaps in reference-guided scaffolded assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_gapfiller (5-align.smk)
    Uses scrubbed reads (not bamtofastq reads) for gap filling.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_GAPFILLER {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/gap2seq:3.1.1a--py311hc84137b_6"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/gap2seq", mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample_id), path(scaffold), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: filled

    script:
    """
    set +e
    Gap2Seq \\
        --scaffolds ${scaffold} \\
        --filled ${sample_id}.fasta \\
        --reads ${r1},${r2} \\
        --threads ${task.cpus} > gap2seq.log 2>&1
    exitcode=\$?
    set -e

    if [ \$exitcode -ne 0 ]; then
        echo "Gap2Seq failed, copying scaffold as filled sequence"
        cp ${scaffold} ${sample_id}.fasta
    fi
    """
}
