/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SELECT_BEST_ASSEMBLY - Choose best assembly between denovo and ref-guided
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: select_best_assembly (6-select-best.smk)
    In auto mode, scores both assemblies and picks the winner.
    In single mode, just copies the appropriate assembly.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SELECT_BEST_ASSEMBLY {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/${sample_id}/best",                     mode: 'copy', pattern: 'assembly.fasta'
    publishDir "${params.outdir}/reports/assembly_selection",             mode: 'copy', pattern: '*.txt'

    input:
    tuple val(sample_id), path(denovo_assembly), path(refguided_assembly), path(denovo_quast), path(refguided_quast)

    output:
    tuple val(sample_id), path("assembly.fasta"),       emit: assembly
    path "${sample_id}.txt",                             emit: report
    path "${sample_id}.method.txt",                      emit: method

    script:
    """
    #!/usr/bin/env python3
    import csv
    import os
    import shutil

    def parse_quast(tsv_path):
        \"\"\"Parse QUAST report.tsv into a dict.\"\"\"
        metrics = {}
        with open(tsv_path) as f:
            reader = csv.reader(f, delimiter='\\t')
            for row in reader:
                if len(row) >= 2:
                    metrics[row[0]] = row[1]
        return metrics

    def score_assembly(stats, expected_length=1040000):
        \"\"\"Calculate quality score (0-100 scale) matching Snakemake select_best_assembly.py.\"\"\"
        score = 0

        n50 = int(stats.get('N50', '0').replace(',', ''))
        num_contigs = int(stats.get('# contigs', '9999').replace(',', ''))
        total_length = int(stats.get('Total length', '0').replace(',', ''))
        gc_percent = float(stats.get('GC (%)', '0'))

        # N50 Score (40 points)
        n50_ratio = n50 / expected_length if expected_length > 0 else 0
        if n50_ratio >= 0.95:
            score += 40
        elif n50_ratio >= 0.5:
            score += 40 * n50_ratio
        else:
            score += 40 * (n50_ratio ** 2)

        # Contig Count Score (30 points)
        if num_contigs <= 2:
            score += 30
        elif num_contigs <= 5:
            score += 25
        elif num_contigs <= 10:
            score += 20
        elif num_contigs <= 20:
            score += 10
        else:
            score += max(0, 10 - (num_contigs - 20))

        # Length Accuracy Score (20 points)
        length_diff = abs(total_length - expected_length)
        if length_diff <= 10000:
            score += 20
        elif length_diff <= 50000:
            score += 20 * (1 - length_diff / 50000)
        else:
            score += max(0, 10 - (length_diff - 50000) / 10000)

        # GC Content Check (10 points)
        if 40.0 <= gc_percent <= 42.0:
            score += 10
        elif 39.0 <= gc_percent <= 43.0:
            score += 8
        else:
            score += max(0, 8 - abs(gc_percent - 41.0))

        return round(score, 2)

    denovo_stats = parse_quast("${denovo_quast}")
    refguided_stats = parse_quast("${refguided_quast}")

    denovo_score = score_assembly(denovo_stats)
    refguided_score = score_assembly(refguided_stats)

    # Minimum margin required for ref-denovo to win
    # De novo is less biased, so prefer it in ties/close calls
    MIN_MARGIN_FOR_REFGUIDED = 2.0

    if refguided_score >= denovo_score + MIN_MARGIN_FOR_REFGUIDED:
        method = "reference-denovo"
        selection_reason = f"ref-denovo wins by {refguided_score - denovo_score:.2f} points (>= {MIN_MARGIN_FOR_REFGUIDED} margin)"
    else:
        method = "denovo"
        if denovo_score > refguided_score:
            selection_reason = f"denovo wins by {denovo_score - refguided_score:.2f} points"
        elif denovo_score == refguided_score:
            selection_reason = f"tie ({denovo_score:.2f} each) - denovo preferred (less biased)"
        else:
            selection_reason = f"ref-denovo ahead by {refguided_score - denovo_score:.2f} points but < {MIN_MARGIN_FOR_REFGUIDED} margin - denovo preferred"

    with open("${sample_id}.txt", 'w') as report:
        report.write(f"Assembly Selection Report\\n")
        report.write(f"========================\\n\\n")
        report.write(f"SELECTED: {method} (Score: {denovo_score if method == 'denovo' else refguided_score:.2f})\\n")
        report.write(f"Reason: {selection_reason}\\n\\n")
        report.write(f"De Novo Assembly:\\n")
        report.write(f"  Score: {denovo_score:.2f}\\n")
        report.write(f"  N50: {denovo_stats.get('N50', 'N/A')} bp\\n")
        report.write(f"  Contigs: {denovo_stats.get('# contigs', 'N/A')}\\n")
        report.write(f"  Length: {denovo_stats.get('Total length', 'N/A')} bp\\n")
        report.write(f"  GC%: {denovo_stats.get('GC (%)', 'N/A')}\\n\\n")
        report.write(f"Reference-Guided Assembly:\\n")
        report.write(f"  Score: {refguided_score:.2f}\\n")
        report.write(f"  N50: {refguided_stats.get('N50', 'N/A')} bp\\n")
        report.write(f"  Contigs: {refguided_stats.get('# contigs', 'N/A')}\\n")
        report.write(f"  Length: {refguided_stats.get('Total length', 'N/A')} bp\\n")
        report.write(f"  GC%: {refguided_stats.get('GC (%)', 'N/A')}\\n\\n")
        report.write(f"Score Difference: {abs(denovo_score - refguided_score):.2f} points\\n")
        report.write(f"Minimum margin for ref-denovo: {MIN_MARGIN_FOR_REFGUIDED} points\\n")

    with open("${sample_id}.method.txt", 'w') as f:
        f.write(method + "\\n")

    if method == "denovo":
        shutil.copy("${denovo_assembly}", "assembly.fasta")
    else:
        shutil.copy("${refguided_assembly}", "assembly.fasta")
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COPY_BEST_ASSEMBLY - Copy single-mode assembly as best
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Used when mode is 'denovo' or 'reference-denovo' (no comparison needed).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COPY_BEST_ASSEMBLY {
    tag "${sample_id}"
    label 'process_low'
    
    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/${sample_id}/best",                     mode: 'copy', pattern: 'assembly.fasta'
    publishDir "${params.outdir}/reports/assembly_selection",             mode: 'copy', pattern: '*.method.txt'

    input:
    tuple val(sample_id), path(assembly)
    val method_name

    output:
    tuple val(sample_id), path("assembly.fasta"), emit: assembly
    path "${sample_id}.method.txt",                emit: method

    script:
    """
    cp ${assembly} assembly.fasta
    echo "${method_name}" > ${sample_id}.method.txt
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ASSEMBLY_METHOD_SUMMARY - Summarize assembly method selections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: assembly_method_summary (6-select-best.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ASSEMBLY_METHOD_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('methods/*')
    val mode

    output:
    path "assembly_methods_summary.txt", emit: summary

    script:
    """
    total=\$(ls methods/*.method.txt 2>/dev/null | wc -l | tr -d ' ')
    denovo_count=\$(grep -rl "^denovo\$" methods/ 2>/dev/null | wc -l | tr -d ' ')
    refguided_count=\$(grep -rl "^reference-denovo\$" methods/ 2>/dev/null | wc -l | tr -d ' ')

    cat > assembly_methods_summary.txt << ENDREPORT
Assembly Method Selection Summary
=================================

ENDREPORT

    if [ "${mode}" = "auto" ]; then
        echo "Mode: auto (automatic selection)" >> assembly_methods_summary.txt
        echo "" >> assembly_methods_summary.txt
        echo "Total samples: \$total" >> assembly_methods_summary.txt
        echo "De novo selected: \$denovo_count samples" >> assembly_methods_summary.txt
        echo "Reference-guided selected: \$refguided_count samples" >> assembly_methods_summary.txt
    else
        echo "Mode: ${mode} (single method)" >> assembly_methods_summary.txt
        echo "All \$total samples used ${mode}" >> assembly_methods_summary.txt
    fi
    """
}
