/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUAST - Assembly quality statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: quast (4-assemble.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process QUAST {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"

    publishDir "${params.outdir}/${sample_id}/denovo/assembly_statistics", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.quast.report.tsv"),   emit: report
    path "${sample_id}.quast.report.html",                          emit: html, optional: true
    path "quast.log",                                             emit: log
    path "*.txt",                                                    emit: txt, optional: true

    script:
    """
    # Preserve QUAST's exit status when logging through tee.
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

    # Move outputs to working dir with sample-prefixed names
    cp quast_out/report.tsv ${sample_id}.quast.report.tsv
    cp quast_out/report.html ${sample_id}.quast.report.html 2>/dev/null || true
    cp quast_out/*.txt . 2>/dev/null || true
    """
}
