/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_COLLATE_COVERAGE - Collate ref-guided coverage TSVs across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_collate_coverage (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_COLLATE_COVERAGE {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}", mode: 'copy', pattern: 'ref-denovo.coverage.tsv'

    input:
    path ('coverages/*')

    output:
    path "ref-denovo.coverage.tsv", emit: coverage

    script:
    """
    first_cov=\$(ls coverages/* | head -1)
    head -1 "\$first_cov" > ref-denovo.coverage.tsv
    for f in coverages/*; do
        tail -n +2 "\$f" >> ref-denovo.coverage.tsv
    done
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_COLLATE_FLAGSTAT - Collate ref-guided alignment stats across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_collate_flagstat (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_COLLATE_FLAGSTAT {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('flagstats/*')

    output:
    path "ref-denovo.alignment_summary.txt", emit: summary

    script:
    """
    echo "Reference-Guided Alignment Summary" > ref-denovo.alignment_summary.txt
    echo "Generated: \$(date)" >> ref-denovo.alignment_summary.txt
    echo "================================" >> ref-denovo.alignment_summary.txt
    echo "" >> ref-denovo.alignment_summary.txt

    for flagstat in flagstats/*; do
        sample=\$(basename "\$flagstat" .flagstat.txt)
        echo "Sample: \$sample" >> ref-denovo.alignment_summary.txt
        cat "\$flagstat" >> ref-denovo.alignment_summary.txt
        echo "" >> ref-denovo.alignment_summary.txt
    done
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REF_SCAFFOLD_SUMMARY - Summarize ref-guided scaffolding results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: ref_scaffold_summary (5-align.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REF_SCAFFOLD_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('results/*')

    output:
    path "ref_scaffold_summary.txt", emit: report

    script:
    """
    echo "RagTag Scaffolding Summary Report (Reference-Guided Mode)" > ref_scaffold_summary.txt
    echo "Generated: \$(date)" >> ref_scaffold_summary.txt
    echo "================================" >> ref_scaffold_summary.txt
    echo "" >> ref_scaffold_summary.txt

    success=\$(grep -l "SCAFFOLD_SUCCESS" results/* 2>/dev/null | wc -l)
    failed=\$(grep -l "SCAFFOLD_FAILED" results/* 2>/dev/null | wc -l)

    echo "Total samples processed: \$((\$success + \$failed))" >> ref_scaffold_summary.txt
    echo "Successfully scaffolded: \$success" >> ref_scaffold_summary.txt
    echo "Failed scaffolding: \$failed" >> ref_scaffold_summary.txt
    echo "" >> ref_scaffold_summary.txt

    if [ \$failed -gt 0 ]; then
        echo "Samples that FAILED scaffolding (using original contigs):" >> ref_scaffold_summary.txt
        grep -l "SCAFFOLD_FAILED" results/* 2>/dev/null | \\
            xargs -I% basename % | \\
            sed 's/^/  - /' >> ref_scaffold_summary.txt
        echo "" >> ref_scaffold_summary.txt
    fi

    if [ \$success -gt 0 ]; then
        echo "Samples that SUCCEEDED scaffolding:" >> ref_scaffold_summary.txt
        grep -l "SCAFFOLD_SUCCESS" results/* 2>/dev/null | \\
            xargs -I% basename % | \\
            sed 's/^/  - /' >> ref_scaffold_summary.txt
    fi
    """
}
