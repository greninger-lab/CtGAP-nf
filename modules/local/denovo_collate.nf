/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DENOVO_COLLATE_COVERAGE - Collate coverage TSVs across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: denovo_collate_coverage (4-assemble.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DENOVO_COLLATE_COVERAGE {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}", mode: 'copy', pattern: 'denovo.coverage.tsv'

    input:
    path ('coverages/*')

    output:
    path "denovo.coverage.tsv", emit: coverage

    script:
    """
    first_cov=\$(ls coverages/* | head -1)
    head -1 "\$first_cov" > denovo.coverage.tsv
    for f in coverages/*; do
        tail -n +2 "\$f" >> denovo.coverage.tsv
    done
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DENOVO_COLLATE_FLAGSTAT - Collate alignment statistics across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: denovo_collate_flagstat (4-assemble.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process DENOVO_COLLATE_FLAGSTAT {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('flagstats/*')

    output:
    path "denovo.alignment_summary.txt", emit: summary

    script:
    """
    echo "De Novo Alignment Summary (vs Reference Set)" > denovo.alignment_summary.txt
    echo "Generated: \$(date)" >> denovo.alignment_summary.txt
    echo "================================" >> denovo.alignment_summary.txt
    echo "" >> denovo.alignment_summary.txt

    for flagstat in flagstats/*; do
        sample=\$(basename "\$flagstat" .flagstat.txt)
        echo "Sample: \$sample" >> denovo.alignment_summary.txt
        cat "\$flagstat" >> denovo.alignment_summary.txt
        echo "" >> denovo.alignment_summary.txt
    done
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCAFFOLD_SUMMARY - Summarize scaffolding results across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: scaffold_summary (4-assemble.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SCAFFOLD_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('results/*')

    output:
    path "scaffold_summary.txt", emit: report

    script:
    """
    echo "RagTag Scaffolding Summary Report" > scaffold_summary.txt
    echo "Generated: \$(date)" >> scaffold_summary.txt
    echo "================================" >> scaffold_summary.txt
    echo "" >> scaffold_summary.txt

    success=\$(grep -l "SCAFFOLD_SUCCESS" results/* 2>/dev/null | wc -l)
    failed=\$(grep -l "SCAFFOLD_FAILED" results/* 2>/dev/null | wc -l)

    echo "Total samples processed: \$((\$success + \$failed))" >> scaffold_summary.txt
    echo "Successfully scaffolded: \$success" >> scaffold_summary.txt
    echo "Failed scaffolding: \$failed" >> scaffold_summary.txt
    echo "" >> scaffold_summary.txt

    if [ \$failed -gt 0 ]; then
        echo "Samples that FAILED scaffolding (using original contigs):" >> scaffold_summary.txt
        grep -l "SCAFFOLD_FAILED" results/* 2>/dev/null | \\
            xargs -I% basename % _result.txt | \\
            sed 's/^/  - /' >> scaffold_summary.txt
        echo "" >> scaffold_summary.txt
    fi

    if [ \$success -gt 0 ]; then
        echo "Samples that SUCCEEDED scaffolding:" >> scaffold_summary.txt
        grep -l "SCAFFOLD_SUCCESS" results/* 2>/dev/null | \\
            xargs -I% basename % _result.txt | \\
            sed 's/^/  - /' >> scaffold_summary.txt
    fi
    """
}
