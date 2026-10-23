/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTER_ASSEMBLY_CONTIGS - Remove non-Ct contigs via minimap2 alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules: filter_assembly_contigs + apply_filtered_assembly
    (6b-filter-contigs.smk)

    Aligns contigs to Ct reference, keeps those with sufficient identity/coverage.
    Outputs the filtered assembly (replacing original) plus a backup and report.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * Step 1: Index assembly with samtools faidx
 */
process FILTER_INDEX_ASSEMBLY {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"    

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path(assembly), path("${assembly}.fai"), emit: indexed

    script:
    """
    samtools faidx ${assembly}
    """
}

/*
 * Step 2: Align contigs to Ct reference with minimap2
 */
process FILTER_MINIMAP2_ALIGN {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/minimap2:2.31--h118bc1c_0"

    input:
    tuple val(sample_id), path(assembly), path(fai)

    output:
    tuple val(sample_id), path(assembly), path(fai), path("assembly.paf"), emit: aligned

    script:
    """
    minimap2 -x asm5 -t ${task.cpus} ${params.ref_genome} ${assembly} > assembly.paf 2> minimap2.log
    """
}

/*
 * Step 3: Filter contigs based on PAF alignment and generate report
 */
process FILTER_APPLY_CONTIGS {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/seqtk:1.5--h577a1d6_1"


    publishDir "${params.outdir}/${sample_id}/best", mode: 'copy', pattern: '*.assembly.*.fasta'
    publishDir "${params.outdir}/${sample_id}/best", mode: 'copy', pattern: '*.filter_report.txt'

    input:
    tuple val(sample_id), path(assembly), path(fai), path(paf)

    output:
    tuple val(sample_id), path("${sample_id}.assembly.filtered.fasta"),  emit: assembly
    path "${sample_id}.assembly.unfiltered.fasta",                        emit: backup
    path "${sample_id}.filter_report.txt",                   emit: report

    script:
    def min_identity = params.filter_min_identity
    def min_coverage = params.filter_min_coverage
    """
    # Parse PAF and identify Ct contigs
    awk -v min_id=${min_identity} -v min_cov=${min_coverage} '
    BEGIN { OFS="\\t" }
    {
        qname = \$1
        qlen = \$2
        matches = \$10
        alnlen = \$11

        identity = (alnlen > 0) ? (matches / alnlen) * 100 : 0
        coverage = (qlen > 0) ? (alnlen / qlen) * 100 : 0

        if (!(qname in best_cov) || coverage > best_cov[qname]) {
            best_cov[qname] = coverage
            best_id[qname] = identity
            best_len[qname] = qlen
        }
    }
    END {
        for (contig in best_cov) {
            if (best_id[contig] >= min_id && best_cov[contig] >= min_cov) {
                print contig, best_len[contig], best_id[contig], best_cov[contig], "KEEP"
            } else {
                print contig, best_len[contig], best_id[contig], best_cov[contig], "REMOVE"
            }
        }
    }' ${paf} > filter_decisions.tsv

    # Check for contigs with NO alignment
    cut -f1 ${fai} | while read contig; do
        if ! grep -q "^\$contig" filter_decisions.tsv; then
            len=\$(grep "^\$contig" ${fai} | cut -f2)
            echo -e "\$contig\\t\$len\\t0\\t0\\tREMOVE_NO_ALIGN"
        fi
    done >> filter_decisions.tsv

    # Extract contigs to keep
    grep "KEEP" filter_decisions.tsv | cut -f1 > keep_contigs.txt || true

    # Create filtered assembly
    if [ -s keep_contigs.txt ]; then
        seqtk subseq ${assembly} keep_contigs.txt > ${sample_id}.assembly.filtered.fasta
    else
        cp ${assembly} ${sample_id}.assembly.filtered.fasta
        echo "WARNING: No contigs passed filtering - keeping original"
    fi

    # Backup original
    cp ${assembly} ${sample_id}.assembly.unfiltered.fasta

    # Generate report
    total_contigs=\$(wc -l < ${fai})
    kept_contigs=\$(grep -c "KEEP" filter_decisions.tsv) || kept_contigs=0
    removed_contigs=\$((total_contigs - kept_contigs))

    original_size=\$(awk '{sum+=\$2} END {print sum}' ${fai})
    filtered_size=\$(seqtk comp ${sample_id}.assembly.filtered.fasta | awk '{sum+=\$2} END {print sum}') || filtered_size=0
    removed_bp=\$((original_size - filtered_size))

    cat > ${sample_id}.filter_report.txt << EOF
ASSEMBLY FILTERING REPORT
=========================
Sample: ${sample_id}
Reference: ${params.ref_genome}
Min identity: ${min_identity}%
Min coverage: ${min_coverage}%

RESULTS
-------
Total contigs: \$total_contigs
Kept contigs: \$kept_contigs
Removed contigs: \$removed_contigs

Original size: \$original_size bp
Filtered size: \$filtered_size bp
Removed: \$removed_bp bp

CONTIG DETAILS
--------------
EOF

    echo -e "Contig\\tLength\\tIdentity\\tCoverage\\tStatus" >> ${sample_id}.filter_report.txt
    cat filter_decisions.tsv >> ${sample_id}.filter_report.txt

    # Cleanup
    rm -f filter_decisions.tsv keep_contigs.txt
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTER_SUMMARY - Summarize filtering across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: filter_summary (6b-filter-contigs.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FILTER_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('reports/*')

    output:
    path "assembly_filter_summary.tsv", emit: summary

    script:
    """
    echo -e "sample\\toriginal_bp\\tfiltered_bp\\tremoved_bp\\ttotal_contigs\\tkept_contigs" > assembly_filter_summary.tsv

    for report in reports/*; do
        sample=\$(grep "^Sample:" "\$report" | head -1 | awk '{print \$2}')
        [ -z "\$sample" ] && sample=\$(basename "\$report" .txt)
        original=\$(grep "^Original size:" "\$report" | head -1 | awk '{print \$3}')
        filtered=\$(grep "^Filtered size:" "\$report" | head -1 | awk '{print \$3}')
        removed_bp=\$(grep "^Removed:" "\$report" | head -1 | awk '{print \$2}')
        total=\$(grep "^Total contigs:" "\$report" | head -1 | awk '{print \$3}')
        kept=\$(grep "^Kept contigs:" "\$report" | head -1 | awk '{print \$3}')
        echo -e "\${sample:-NA}\\t\${original:-NA}\\t\${filtered:-NA}\\t\${removed_bp:-NA}\\t\${total:-NA}\\t\${kept:-NA}" >> assembly_filter_summary.tsv
    done

    sort -t\$'\\t' -k1,1 -o assembly_filter_summary.tsv <(head -1 assembly_filter_summary.tsv) <(tail -n +2 assembly_filter_summary.tsv | sort -t\$'\\t' -k1,1)
    """
}
