/* Gate assembly on mapped coverage and the actual FASTQs destined for Shovill. */
process CHECK_ASSEMBLY_COVERAGE {
    tag "${sample_id}:${method}"
    label 'process_low'
    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports/coverage_checks", mode: 'copy', pattern: '*.coverage_check.tsv'
    publishDir "${params.outdir}/02_reports/errors", mode: 'copy', pattern: '*.coverage_error.txt'
    publishDir "${params.outdir}/01_per_sample/${sample_id}", mode: 'copy', pattern: '*.coverage_error.txt'

    input:
    tuple val(sample_id), path(r1), path(r2), path(coverage)
    val method
    path coverage_checker

    output:
    tuple val(sample_id), env(PASSES), emit: decision
    path "${sample_id}.${method}.coverage_check.tsv", emit: metrics
    path "${sample_id}.${method}.coverage_error.txt", emit: error_report, optional: true

    script:
    """
    PASSES=\$(python3 '${coverage_checker}' \\
        --sample '${sample_id}' --method '${method}' \\
        --r1 '${r1}' --r2 '${r2}' --coverage '${coverage}' \\
        --minimum-depth ${params.min_assembly_depth} \\
        --metrics '${sample_id}.${method}.coverage_check.tsv' \\
        --error-report '${sample_id}.${method}.coverage_error.txt')
    """
}

process ASSEMBLY_COVERAGE_SUMMARY {
    tag 'all_samples'
    label 'process_low'
    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/02_reports/summaries", mode: 'copy'

    input:
    path ('checks/*')

    output:
    path 'assembly_coverage_checks.tsv', emit: checks
    path 'assembly_errors.tsv', emit: errors
    path 'assembly_coverage_summary.txt', emit: summary

    script:
    """
    #!/usr/bin/env python3
    import csv
    from collections import defaultdict
    from pathlib import Path

    rows = []
    for path in sorted(Path('checks').glob('*.coverage_check.tsv')):
        with path.open() as handle:
            rows.extend(csv.DictReader(handle, delimiter='\\t'))
    fields = list(rows[0])
    for filename, selected in [('assembly_coverage_checks.tsv', rows),
                                ('assembly_errors.tsv', [r for r in rows if r['status'] != 'PASS'])]:
        with open(filename, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\\t')
            writer.writeheader()
            writer.writerows(selected)
    samples = defaultdict(list)
    for row in rows:
        samples[row['sample']].append(row)
    rejected = sorted(sample for sample, checks in samples.items() if not any(r['status'] == 'PASS' for r in checks))
    with open('assembly_coverage_summary.txt', 'w') as handle:
        handle.write('CtGAP Assembly Coverage Summary\\n')
        handle.write('Minimum mapped mean depth: ${params.min_assembly_depth}x\\n')
        handle.write(f'Samples checked: {len(samples)}\\n')
        handle.write(f'Samples with an eligible assembly path: {len(samples) - len(rejected)}\\n')
        handle.write(f'Samples with all requested paths skipped: {len(rejected)}\\n')
        for sample in rejected:
            reasons = ', '.join(r['assembly_path'] + ': ' + r['reason'] for r in samples[sample])
            handle.write(f'{sample}: ASSEMBLY_SKIPPED ({reasons})\\n')
    """
}
