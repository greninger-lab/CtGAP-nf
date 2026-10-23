/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_SCAFFOLD - Scaffold reference-guided assembly with RagTag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_scaffold (5-align.smk)
    Falls back to original contigs if RagTag fails.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_SCAFFOLD {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/ragtag:2.1.0--pyhdfd78af_2"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/scaffold", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("ragtag.scaffold.fasta"),            emit: scaffold
    path "${sample_id}.scaffold_result.txt",                         emit: result

    script:
    """
    set +e
    ragtag.py scaffold \\
        ${params.scaffold_ref} \\
        ${contigs} \\
        -o ragtag_out \\
        -f ${params.ragtag_min_len} \\
        -r \\
        -u \\
        -t ${task.cpus} \\
        2> ragtag.log
    exitcode=\$?
    set -e

    if [ \$exitcode -ne 0 ] || [ ! -f ragtag_out/ragtag.scaffold.fasta ]; then
        echo "WARNING: RagTag scaffolding failed for ${sample_id}. Using original contigs."
        echo "SCAFFOLD_FAILED" > ${sample_id}.scaffold_result.txt
        cp ${contigs} ragtag.scaffold.fasta
    else
        echo "SCAFFOLD_SUCCESS" > ${sample_id}.scaffold_result.txt
        cp ragtag_out/ragtag.scaffold.fasta .
    fi
    """
}
