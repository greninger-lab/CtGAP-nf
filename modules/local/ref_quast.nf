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

    script:
    """
    quast \\
        ${assembly} \\
        -t ${task.cpus} \\
        -o quast_out 2>&1 | tee quast.log

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

    script:
    """
    quast \\
        ${assembly} \\
        -R ${params.scaffold_ref} \\
        -t ${task.cpus} \\
        -o quast_out 2>&1 | tee quast.log

    cp quast_out/report.tsv ${sample_id}.ref_quast_scaffolded.report.tsv
    """
}
