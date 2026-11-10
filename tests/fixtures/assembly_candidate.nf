process CREATE_CANDIDATE {
    input:
    tuple val(sample_id), path(r1), path(r2)
    val method

    output:
    tuple val(sample_id), path("${sample_id}.${method}.fasta"), emit: assembly
    tuple val(sample_id), path("${sample_id}.${method}.quast.tsv"), emit: report

    script:
    """
    printf '>contig_%s\\nACGT\\n' '${method}' > '${sample_id}.${method}.fasta'
    # Equal, intentionally poor scores exercise the single-path fallback independently of scoring.
    printf 'Assembly\\ttest\\nN50\\t4\\n# contigs\\t9999\\nTotal length\\t4\\nGC (%%)\\t0\\n' > '${sample_id}.${method}.quast.tsv'
    """
}
