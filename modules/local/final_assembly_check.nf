process CHECK_FINAL_ASSEMBLY {
    tag "${sample_id}"
    label 'process_low'
    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports/assembly_checks", mode: 'copy', pattern: '*.final_assembly_check.tsv'
    publishDir "${params.outdir}/02_reports/errors", mode: 'copy', pattern: '*.final_assembly_error.txt'
    publishDir "${params.outdir}/01_per_sample/${sample_id}", mode: 'copy', pattern: '*.final_assembly_error.txt'

    input:
    tuple val(sample_id), path(assembly)
    path checker

    output:
    tuple val(sample_id), env(PASSES), emit: decision
    path "${sample_id}.final_assembly_check.tsv", emit: metrics
    path "${sample_id}.final_assembly_error.txt", emit: error_report, optional: true

    script:
    """
    PASSES=\$(python3 '${checker}' \\
        --sample '${sample_id}' --assembly '${assembly}' \\
        --minimum-length ${params.quast_min_contig} \\
        --metrics '${sample_id}.final_assembly_check.tsv' \\
        --error-report '${sample_id}.final_assembly_error.txt')
    """
}

process FINAL_ASSEMBLY_SUMMARY {
    tag 'all_samples'
    label 'process_low'
    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/02_reports/summaries", mode: 'copy'

    input:
    path ('checks/*')

    output:
    path 'final_assembly_checks.tsv', emit: checks
    path 'final_assembly_errors.tsv', emit: errors

    script:
    """
    #!/usr/bin/env python3
    import csv
    from pathlib import Path

    rows = []
    for path in sorted(Path('checks').glob('*.final_assembly_check.tsv')):
        with path.open() as handle:
            rows.extend(csv.DictReader(handle, delimiter='\\t'))
    for filename, selected in [('final_assembly_checks.tsv', rows),
                                ('final_assembly_errors.tsv', [r for r in rows if r['status'] != 'PASS'])]:
        with open(filename, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\\t')
            writer.writeheader()
            writer.writerows(selected)
    """
}
