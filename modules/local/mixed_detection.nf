/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MIXED STRAIN DETECTION MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules from 11-mixed-detection.smk:
      mixed_extract_main_contig, mixed_index_assembly, mixed_map_reads,
      mixed_call_variants, mixed_detection_summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * Extract largest contig from assembly for mixed strain analysis
 */
process MIXED_EXTRACT_MAIN_CONTIG {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1" 

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("main_contig.fasta"), emit: main_contig

    script:
    """
    # Find the largest contig and extract it using the .fai index (most reliable)
    samtools faidx ${assembly}

    # The .fai file has: name, length, offset, linebases, linewidth
    # Sort by length (column 2) descending and take the first
    largest_contig=\$(sort -t\$'\\t' -k2,2nr ${assembly}.fai | head -1 | cut -f1)

    if [ -z "\$largest_contig" ]; then
        echo "WARNING: Could not determine largest contig" >&2
        cp ${assembly} main_contig.fasta
    else
        samtools faidx ${assembly} "\$largest_contig" > main_contig.fasta
    fi
    """
}

/*
 * Step 1: Build bowtie2 index on main contig
 */
process MIXED_BOWTIE_INDEX {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1"

    input:
    tuple val(sample_id), path(main_contig), path(r1), path(r2)

    output:
    tuple val(sample_id), path(main_contig), path(r1), path(r2), path("assembly.*.bt2"), emit: indexed

    script:
    """
    bowtie2-build ${main_contig} assembly > /dev/null 2>&1
    """
}

/*
 * Step 2: Align reads to main contig with bowtie2
 */
process MIXED_BOWTIE_ALIGN {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1"

    input:
    tuple val(sample_id), path(main_contig), path(r1), path(r2), path(bt2_idx)

    output:
    tuple val(sample_id), path(main_contig), path("${sample_id}.mixed.sam"), emit: aligned

    script:
    """
    bowtie2 \\
        -x assembly \\
        -1 ${r1} \\
        -2 ${r2} \\
        --threads ${task.cpus} \\
        -S ${sample_id}.mixed.sam \\
        2> bowtie2.log
    """
}

/*
 * Step 3: Sort, index BAM and index FASTA with samtools
 */
process MIXED_SAMTOOLS_SORT {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

    input:
    tuple val(sample_id), path(main_contig), path(sam)

    output:
    tuple val(sample_id), path(main_contig), path("${main_contig}.fai"), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted

    script:
    """
    samtools view -bSh ${sam} | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam

    samtools faidx ${main_contig}
    """
}

/*
 * Step 4: Call variants with bcftools and extract het sites
 */
process MIXED_BCFTOOLS_CALL {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/bcftools:1.24--h118bc1c_4"

    publishDir "${params.outdir}/${sample_id}/mixed_detection", mode: 'copy'

    input:
    tuple val(sample_id), path(main_contig), path(fai), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.variants.vcf"),       emit: vcf
    tuple val(sample_id), path("${sample_id}.het_sites.tsv"),      emit: het_sites
    tuple val(sample_id), path("${sample_id}.mixed_summary.txt"), emit: summary

    script:
    def min_af = params.mixed_min_af
    def max_af = params.mixed_max_af
    def min_depth = params.mixed_min_depth
    def min_base_qual = params.mixed_min_base_qual
    def min_strand_reads = params.mixed_min_strand_reads
    def het_threshold = params.mixed_het_threshold
    """
    # Call variants
    bcftools mpileup \\
        -f ${main_contig} \\
        -a AD,ADF,ADR,DP \\
        -Q ${min_base_qual} \\
        --max-depth 1000 \\
        ${bam} 2> mpileup.log | \\
    bcftools call \\
        -mv \\
        -Ov \\
        -o ${sample_id}.variants.vcf 2> call.log

    # Extract heterozygous-like sites with strand bias filter
    echo -e "contig\\tposition\\tref\\talt\\tdepth\\tref_count\\talt_count\\talt_freq\\tfwd_alt\\trev_alt\\tstrand_bias" > ${sample_id}.het_sites.tsv

    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%DP\\t[%AD]\\t[%ADF]\\t[%ADR]\\n' \\
        ${sample_id}.variants.vcf 2> /dev/null | \\
    awk -F'\\t' -v min_af=${min_af} -v max_af=${max_af} -v min_dp=${min_depth} -v min_strand=${min_strand_reads} '
    BEGIN { OFS="\\t" }
    {
        split(\$6, ad, ",")
        ref_count = ad[1]
        alt_count = ad[2]
        depth = \$5

        split(\$7, adf, ",")
        split(\$8, adr, ",")
        fwd_alt = (length(adf) > 1) ? adf[2] : 0
        rev_alt = (length(adr) > 1) ? adr[2] : 0

        if (depth >= min_dp && alt_count > 0) {
            af = alt_count / depth
            if (af >= min_af && af <= max_af) {
                if (fwd_alt >= min_strand && rev_alt >= min_strand) {
                    strand_bias = "PASS"
                    print \$1, \$2, \$3, \$4, depth, ref_count, alt_count, af, fwd_alt, rev_alt, strand_bias
                }
            }
        }
    }' >> ${sample_id}.het_sites.tsv

    # Count het sites
    het_count=\$(tail -n +2 ${sample_id}.het_sites.tsv | wc -l | tr -d ' ')

    # Mean AF
    mean_af=\$(tail -n +2 ${sample_id}.het_sites.tsv | awk '{ sum += \$8; n++ } END { if (n > 0) printf "%.3f", sum/n; else print "N/A" }')

    # Main contig info
    main_contig_name=\$(head -1 ${main_contig} | sed 's/>//' | cut -d' ' -f1)
    main_contig_size=\$(tail -n +2 ${main_contig} | tr -d '\\n' | wc -c | tr -d ' ')

    # Position range
    if [ "\$het_count" -gt 0 ]; then
        position_range=\$(tail -n +2 ${sample_id}.het_sites.tsv | awk '
            BEGIN { min=999999999; max=0 }
            { if (\$2 < min) min=\$2; if (\$2 > max) max=\$2 }
            END { if (max > 0) print max - min; else print 0 }
        ')
    else
        position_range=0
    fi

    # Generate summary
    cat > ${sample_id}.mixed_summary.txt << EOF
MIXED STRAIN DETECTION SUMMARY
==============================
Sample: ${sample_id}
Analysis: Variant heterozygosity on main chromosome (largest contig)

PARAMETERS
----------
Min allele frequency: ${min_af}
Max allele frequency: ${max_af}
Min read depth: ${min_depth}
Min base quality: ${min_base_qual}
Min reads per strand: ${min_strand_reads}
Het site threshold: ${het_threshold}

RESULTS
-------
Main contig analyzed: \$main_contig_name (\$main_contig_size bp)
Heterozygous-like sites: \$het_count
Mean allele frequency: \$mean_af
Position range: \$position_range bp

INTERPRETATION
--------------
EOF

    if [ "\$het_count" -gt ${het_threshold} ]; then
        echo "STATUS: MIXED_STRAIN_WARNING" >> ${sample_id}.mixed_summary.txt
        echo "" >> ${sample_id}.mixed_summary.txt
        echo "HIGH HETEROZYGOSITY DETECTED" >> ${sample_id}.mixed_summary.txt
        echo "Found \$het_count sites with intermediate allele frequencies." >> ${sample_id}.mixed_summary.txt
        if [ "\$position_range" -gt 100000 ]; then
            echo "Distribution: GENOME-WIDE (variants span \$position_range bp)" >> ${sample_id}.mixed_summary.txt
        else
            echo "Distribution: LOCALIZED (variants within \$position_range bp region)" >> ${sample_id}.mixed_summary.txt
        fi
    elif [ "\$het_count" -gt 10 ]; then
        echo "STATUS: LOW_HETEROZYGOSITY" >> ${sample_id}.mixed_summary.txt
        echo "" >> ${sample_id}.mixed_summary.txt
        echo "Moderate heterozygosity detected (\$het_count sites)." >> ${sample_id}.mixed_summary.txt
    else
        echo "STATUS: SINGLE_STRAIN" >> ${sample_id}.mixed_summary.txt
        echo "" >> ${sample_id}.mixed_summary.txt
        echo "Low heterozygosity (\$het_count sites) consistent with single strain." >> ${sample_id}.mixed_summary.txt
    fi
    """
}

/*
 * Collate mixed detection results across all samples
 */
process MIXED_DETECTION_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('summaries/*')
    path ('het_sites/*')

    output:
    path "mixed_detection_summary.tsv", emit: summary_tsv
    path "mixed_detection_report.txt",  emit: report

    script:
    """
    # Build TSV summary
    echo -e "sample\\thet_sites\\tmean_af\\tstatus" > mixed_detection_summary.tsv

    for summary in summaries/*.mixed_summary.txt; do
        sample=\$(basename "\$summary" .mixed_summary.txt)
        het_sites=\$(grep "Heterozygous-like sites:" "\$summary" | awk '{print \$NF}')
        het_sites=\${het_sites:-0}
        mean_af=\$(grep "Mean allele frequency:" "\$summary" | awk '{print \$NF}')
        mean_af=\${mean_af:-N/A}

        if grep -q "STATUS: MIXED_STRAIN_WARNING" "\$summary"; then
            status="MIXED_WARNING"
        elif grep -q "STATUS: LOW_HETEROZYGOSITY" "\$summary"; then
            status="LOW_HET"
        elif grep -q "STATUS: SINGLE_STRAIN" "\$summary"; then
            status="SINGLE"
        else
            status="Unknown"
        fi

        echo -e "\${sample}\\t\${het_sites}\\t\${mean_af}\\t\${status}" >> mixed_detection_summary.tsv
    done

    # Sort by sample name (skip header)
    head -1 mixed_detection_summary.tsv > tmp_header.tsv
    tail -n +2 mixed_detection_summary.tsv | sort -t\$'\\t' -k1,1 > tmp_body.tsv
    cat tmp_header.tsv tmp_body.tsv > mixed_detection_summary.tsv
    rm -f tmp_header.tsv tmp_body.tsv

    # Build report
    total=\$(tail -n +2 mixed_detection_summary.tsv | wc -l | tr -d ' ')
    mixed_count=\$(tail -n +2 mixed_detection_summary.tsv | awk -F'\\t' '\$4=="MIXED_WARNING"' | wc -l | tr -d ' ')
    low_het_count=\$(tail -n +2 mixed_detection_summary.tsv | awk -F'\\t' '\$4=="LOW_HET"' | wc -l | tr -d ' ')
    single_count=\$(tail -n +2 mixed_detection_summary.tsv | awk -F'\\t' '\$4=="SINGLE"' | wc -l | tr -d ' ')

    cat > mixed_detection_report.txt << ENDREPORT
======================================================================
MIXED STRAIN DETECTION REPORT (CtGAP-nf)
======================================================================
Generated: \$(date '+%Y-%m-%d %H:%M:%S')
Total samples analyzed: \$total

----------------------------------------
SUMMARY
----------------------------------------
Single strain (low het):     \$single_count samples
Moderate heterozygosity:     \$low_het_count samples
Mixed strain warning:        \$mixed_count samples

ENDREPORT

    if [ "\$mixed_count" -gt 0 ]; then
        echo "----------------------------------------" >> mixed_detection_report.txt
        echo "SAMPLES WITH MIXED STRAIN WARNING" >> mixed_detection_report.txt
        echo "----------------------------------------" >> mixed_detection_report.txt
        tail -n +2 mixed_detection_summary.tsv | awk -F'\\t' '\$4=="MIXED_WARNING" {printf "  %s: %s het sites (mean AF: %s)\\n", \$1, \$2, \$3}' >> mixed_detection_report.txt
        echo "" >> mixed_detection_report.txt
    fi

    echo "----------------------------------------" >> mixed_detection_report.txt
    echo "ALL SAMPLES" >> mixed_detection_report.txt
    echo "----------------------------------------" >> mixed_detection_report.txt
    printf "%-30s %10s %10s %15s\\n" "Sample" "Het Sites" "Mean AF" "Status" >> mixed_detection_report.txt
    echo "----------------------------------------------------------------------" >> mixed_detection_report.txt
    tail -n +2 mixed_detection_summary.tsv | sort -t\$'\\t' -k2,2nr | awk -F'\\t' '{printf "%-30s %10s %10s %15s\\n", \$1, \$2, \$3, \$4}' >> mixed_detection_report.txt
    echo "" >> mixed_detection_report.txt
    echo "======================================================================" >> mixed_detection_report.txt
    """
}
