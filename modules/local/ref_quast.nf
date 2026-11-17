/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_QUAST - Assembly quality statistics (reference-free) for ref-guided assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_quast (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_QUAST {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.ref_quast.report.tsv"), emit: report
    path "quast.log", emit: log

    script:
    """
    set -o pipefail

    quast \\
        ${assembly} \\
        --min-contig ${params.quast_min_contig} \\
        -t ${task.cpus} \\
        -o quast_out 2>&1 | tee quast.log

    if [ ! -s quast_out/report.tsv ]; then
        printf '%s\\n' '${task.process}: no non-empty report.tsv was produced for ${sample_id}; inspect quast.log and the assembly contig lengths.' >&2
        exit 1
    fi

    cp quast_out/report.tsv ${sample_id}.ref_quast.report.tsv
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_QUAST_SCAFFOLDED - Assembly stats vs scaffold reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_quast_scaffolded (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_QUAST_SCAFFOLDED {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"

    publishDir "${params.outdir}/${sample_id}/ref-denovo/assembly_statistics_vs_scaffold", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.ref_quast_scaffolded.report.tsv"), emit: report
    path "quast.log", emit: log

    script:
    """
    set -o pipefail

    quast \\
        ${assembly} \\
        -R ${params.scaffold_ref} \\
        --min-contig ${params.quast_min_contig} \\
        -t ${task.cpus} \\
        -o quast_out 2>&1 | tee quast.log

    if [ ! -s quast_out/report.tsv ]; then
        printf '%s\\n' '${task.process}: no non-empty report.tsv was produced for ${sample_id}; inspect quast.log and the assembly contig lengths.' >&2
        exit 1
    fi

    cp quast_out/report.tsv ${sample_id}.ref_quast_scaffolded.report.tsv
    """
}
