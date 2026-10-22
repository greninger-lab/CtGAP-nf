/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAKTA_ANNOTATE - Genome annotation using Bakta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules: bakta_annotate (10-annotate.smk)
    Only runs for samples that pass coverage threshold.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CHECK_ANNOTATION_COVERAGE {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    input:
    tuple val(sample_id), path(quast_report)

    output:
    tuple val(sample_id), env(PASSES), emit: result

    script:
    def min_coverage = params.annotation_min_coverage
    def ref_size = params.ct_reference_size
    """
    # Extract N50 from QUAST report (exact match to avoid N50/N75 confusion)
    n50=\$(awk -F'\\t' '\$1 == "N50" {gsub(",","",\$2); print \$2; exit}' ${quast_report})
    n50=\${n50:-0}

    # Calculate coverage percentage
    coverage_pct=\$(awk "BEGIN {printf \\"%.1f\\", (\$n50 / ${ref_size}) * 100}")

    if awk "BEGIN {exit !(\$n50 / ${ref_size} * 100 >= ${min_coverage})}"; then
        PASSES="true"
    else
        PASSES="false"
    fi
    """
}

process BAKTA_ANNOTATE {
    tag "${sample_id}"
    label 'process_high'

    container "quay.io/biocontainers/bakta:1.12.2--pyhdfd78af_0"

    publishDir "${params.outdir}/${sample_id}/annotation", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.gff3"),  emit: gff
    tuple val(sample_id), path("${sample_id}.gbff"),  emit: gbk
    tuple val(sample_id), path("${sample_id}.faa"),   emit: faa
    tuple val(sample_id), path("${sample_id}.ffn"),   emit: ffn
    tuple val(sample_id), path("${sample_id}.tsv"),   emit: tsv
    tuple val(sample_id), path("${sample_id}.txt"),   emit: txt

    script:
    """
    bakta \\
        --db ${params.bakta_db} \\
        --output . \\
        --prefix ${sample_id} \\
        --genus ${params.bakta_genus} \\
        --species ${params.bakta_species} \\
        --gram ${params.bakta_gram} \\
        --threads ${task.cpus} \\
        --force \\
        ${assembly}
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATION_SUMMARY - Collate annotation results across all annotated samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: annotation_summary (10-annotate.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANNOTATION_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('txt_files/*')
    path coverage_check

    output:
    path "annotation_summary.tsv",  emit: summary
    path "annotation_report.txt",   emit: report

    script:
    """
    # Helper to extract a value from bakta txt output
    extract_val() {
        local file="\$1" pattern="\$2"
        grep -i "\$pattern" "\$file" | head -1 | grep -oP '\\d[\\d,]*' | head -1 | tr -d ','
    }

    # Build TSV
    echo -e "sample\\tcontigs\\tsize_bp\\tgc_percent\\tcds_count\\ttrna_count\\trrna_count\\tncrna_count\\thypothetical_count\\tpseudogene_count" > annotation_summary.tsv

    sample_count=0
    for txt in txt_files/*.txt; do
        [ ! -f "\$txt" ] && continue
        sample=\$(basename "\$txt" .txt)
        contigs=\$(extract_val "\$txt" "Sequences:")
        size_bp=\$(extract_val "\$txt" "Size:")
        gc_pct=\$(grep -i "GC:" "\$txt" | head -1 | grep -oP '[\\d.]+' | head -1)
        cds=\$(extract_val "\$txt" "CDSs:")
        trnas=\$(extract_val "\$txt" "tRNAs:")
        rrnas=\$(extract_val "\$txt" "rRNAs:")
        ncrnas=\$(extract_val "\$txt" "ncRNAs:")
        hypo=\$(extract_val "\$txt" "hypotheticals:")
        pseudo=\$(extract_val "\$txt" "pseudogenes:")
        echo -e "\${sample}\\t\${contigs:-}\\t\${size_bp:-}\\t\${gc_pct:-}\\t\${cds:-}\\t\${trnas:-}\\t\${rrnas:-}\\t\${ncrnas:-}\\t\${hypo:-}\\t\${pseudo:-}" >> annotation_summary.tsv
        sample_count=\$((sample_count + 1))
    done

    # Sort by sample name
    head -1 annotation_summary.tsv > tmp_ann.tsv
    tail -n +2 annotation_summary.tsv | sort >> tmp_ann.tsv
    mv tmp_ann.tsv annotation_summary.tsv

    # Build report
    cat > annotation_report.txt << ENDHEADER
======================================================================
CtGAP-nf GENOME ANNOTATION SUMMARY (Bakta)
======================================================================
Generated: \$(date '+%Y-%m-%d %H:%M:%S')
Samples annotated: \$sample_count

----------------------------------------
COVERAGE CHECK RESULTS
----------------------------------------
ENDHEADER

    cat ${coverage_check} >> annotation_report.txt
    echo "" >> annotation_report.txt

    if [ "\$sample_count" -gt 0 ]; then
        # CDS stats
        cds_vals=\$(tail -n +2 annotation_summary.tsv | cut -f5 | grep -v '^\$' | sort -n)
        if [ -n "\$cds_vals" ]; then
            cds_min=\$(echo "\$cds_vals" | head -1)
            cds_max=\$(echo "\$cds_vals" | tail -1)
            cds_count=\$(echo "\$cds_vals" | wc -l | tr -d ' ')
            cds_median=\$(echo "\$cds_vals" | sed -n "\$((cds_count / 2 + 1))p")
            echo "----------------------------------------" >> annotation_report.txt
            echo "ANNOTATION STATISTICS" >> annotation_report.txt
            echo "----------------------------------------" >> annotation_report.txt
            echo "CDS count range: \$cds_min - \$cds_max" >> annotation_report.txt
            echo "Median CDS count: \$cds_median" >> annotation_report.txt
            echo "" >> annotation_report.txt
        fi

        echo "----------------------------------------" >> annotation_report.txt
        echo "PER-SAMPLE DETAILS" >> annotation_report.txt
        echo "----------------------------------------" >> annotation_report.txt
        tail -n +2 annotation_summary.tsv | sort | while IFS=\$'\\t' read -r s contigs size gc cds trnas rrnas ncrnas hypo pseudo; do
            echo "" >> annotation_report.txt
            echo "\${s}:" >> annotation_report.txt
            echo "  CDSs: \${cds:-N/A}" >> annotation_report.txt
            echo "  tRNAs: \${trnas:-N/A}" >> annotation_report.txt
            echo "  rRNAs: \${rrnas:-N/A}" >> annotation_report.txt
            echo "  Hypotheticals: \${hypo:-N/A}" >> annotation_report.txt
        done
    else
        echo "No samples met the coverage threshold for annotation." >> annotation_report.txt
    fi

    echo "" >> annotation_report.txt
    echo "======================================================================" >> annotation_report.txt
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATION_COVERAGE_CHECK - Generate coverage check TSV for all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process ANNOTATION_COVERAGE_CHECK {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path ('quast_reports/*')

    output:
    path "coverage_check.tsv", emit: coverage_check

    script:
    def min_coverage = params.annotation_min_coverage
    def ref_size = params.ct_reference_size
    """
    echo -e "sample\\tn50\\tcoverage_pct\\tpasses_threshold" > coverage_check.tsv

    for report in quast_reports/*; do
        sample=\$(basename "\$report" .tsv)
        sample=\${sample//.report/}
        [ "\$sample" = "report" ] && sample="unknown"

        n50=\$(grep "^N50" "\$report" | cut -f2 | tr -d ',')
        n50=\${n50:-0}

        if [ "\$n50" -gt 0 ]; then
            coverage_pct=\$(awk "BEGIN {printf \\"%.1f\\", (\$n50 / ${ref_size}) * 100}")
        else
            coverage_pct="0.0"
        fi

        passes=\$(awk "BEGIN {if (\$n50 / ${ref_size} * 100 >= ${min_coverage}) print \\"YES\\"; else print \\"NO\\"}")

        echo -e "\${sample}\\t\${n50}\\t\${coverage_pct}\\t\${passes}" >> coverage_check.tsv
    done
    """
}
