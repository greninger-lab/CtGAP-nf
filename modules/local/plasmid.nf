/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PLASMID ASSEMBLY MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules from plasmid_dedicated.smk:
      index_reference_plasmid, map_reads_to_plasmid, extract_plasmid_reads,
      unicycler_plasmid_assembly, analyze_plasmid_assembly, reorient_plasmid,
      blast_plasmid_dedicated, mlst_plasmid_dedicated, collate_plasmid_results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// INDEX_REFERENCE_PLASMID removed — pre-built bt2 index used directly
// from mount point via params.plasmid_bt2_prefix

/*
 * Step 1: Align scrubbed reads to reference plasmid with bowtie2
 */
process PLASMID_BOWTIE_ALIGN {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/bowtie2:2.5.5--h63e9258_1"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid.sam"), emit: sam

    script:
    """
    bowtie2 \\
        -x ${params.plasmid_bt2_prefix} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --threads ${task.cpus} \\
        --sensitive-local \\
        -S ${sample_id}.plasmid.sam \\
        2> bowtie2.log
    """
}

/*
 * Step 2: Filter mapped reads, sort, index, and generate stats with samtools
 */
process PLASMID_SAMTOOLS_SORT {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

    publishDir "${params.outdir}/${sample_id}/plasmid/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid_mapped.bam"), path("${sample_id}.plasmid_mapped.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.mapping_stats.txt"),                                                emit: stats

    script:
    """
    samtools view -bS -F 4 ${sam} | \\
    samtools sort -@ ${task.cpus} -o ${sample_id}.plasmid_mapped.bam -

    samtools index ${sample_id}.plasmid_mapped.bam

    echo "Plasmid Read Mapping Statistics" > ${sample_id}.mapping_stats.txt
    echo "Sample: ${sample_id}" >> ${sample_id}.mapping_stats.txt
    echo "Reference: pCT (D strain)" >> ${sample_id}.mapping_stats.txt
    echo "================================" >> ${sample_id}.mapping_stats.txt
    samtools flagstat ${sample_id}.plasmid_mapped.bam >> ${sample_id}.mapping_stats.txt
    """
}

/*
 * Extract plasmid read pairs from mapped BAM
 */
process EXTRACT_PLASMID_READS {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

    publishDir "${params.outdir}/${sample_id}/plasmid/reads", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid_R1.fastq.gz"), path("${sample_id}.plasmid_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.read_count.txt"),                                                 emit: count

    script:
    """
    # Extract properly paired reads (-f 3)
    samtools view -bf 3 ${bam} | \\
    samtools sort -n -@ ${task.cpus} - | \\
    samtools fastq -@ ${task.cpus} \\
        -1 ${sample_id}.plasmid_R1.fastq.gz \\
        -2 ${sample_id}.plasmid_R2.fastq.gz \\
        -0 /dev/null \\
        -s /dev/null \\
        -n \\
        2> extract.log

    # Count extracted reads
    num_reads=\$(gzip -dc ${sample_id}.plasmid_R1.fastq.gz | wc -l | awk '{print \$1/4}')
    echo "Extracted plasmid read pairs: \$num_reads" > ${sample_id}.read_count.txt
    echo "\$num_reads" >> ${sample_id}.read_count.txt
    """
}

/*
 * Unicycler assembly of plasmid reads
 */
process UNICYCLER_PLASMID {
    tag "${sample_id}"
    label 'process_high'

    container "quay.io/biocontainers/unicycler:0.5.1--py312hdcc493e_5"

    publishDir "${params.outdir}/${sample_id}/plasmid/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), path(count_file)

    output:
    tuple val(sample_id), path("assembly.fasta"), path("assembly.gfa"), emit: assembly

    script:
    def min_reads = params.plasmid_min_reads
    """
    read_count=\$(tail -1 ${count_file})

    if [ "\$read_count" -lt "${min_reads}" ]; then
        echo "Insufficient plasmid reads (\$read_count < ${min_reads})"
        echo ">no_plasmid_detected" > assembly.fasta
        echo "NNNNNNNNNN" >> assembly.fasta
        touch assembly.gfa
    else
        set +e
        unicycler \\
            -1 ${r1} \\
            -2 ${r2} \\
            -o uc_out \\
            -t ${task.cpus} \\
            --min_fasta_length 1000 \\
            --verbosity 2 \\
            > unicycler.log 2>&1
        exitcode=\$?
        set -e

        if [ \$exitcode -ne 0 ] || [ ! -f uc_out/assembly.fasta ]; then
            echo ">assembly_failed" > assembly.fasta
            echo "NNNNNNNNNN" >> assembly.fasta
            touch assembly.gfa
        else
            cp uc_out/assembly.fasta .
            cp uc_out/assembly.gfa .
        fi
    fi
    """
}

/*
 * Analyze plasmid assembly - check size, circularity, select best contig
 */
process ANALYZE_PLASMID {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/${sample_id}/plasmid", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(graph)
    tuple val(sample_id_stats), path(mapping_stats)
    tuple val(sample_id_count), path(read_count)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid.fasta"),        emit: plasmid
    tuple val(sample_id), path("${sample_id}.plasmid_analysis.txt"), emit: analysis

    script:
    def expected = params.plasmid_expected_size
    def min_size = params.plasmid_min_size
    def max_size = params.plasmid_max_size
    """
    #!/usr/bin/env python3
    import os

    sample = "${sample_id}"
    expected_size = ${expected}
    min_size = ${min_size}
    max_size = ${max_size}

    # Parse assembly
    contigs = []
    with open("${assembly}") as f:
        header = ""
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq:
                    contigs.append((header, seq, len(seq)))
                header = line[1:]
                seq = ""
            else:
                seq += line
        if header and seq:
            contigs.append((header, seq, len(seq)))

    # Check for circular contigs in GFA
    circular_contigs = set()
    gfa = "${graph}"
    if os.path.exists(gfa) and os.path.getsize(gfa) > 0:
        with open(gfa) as f:
            for line in f:
                if line.startswith("L") and "circular" in line.lower():
                    circular_contigs.add(line.split()[1])

    NL = chr(10)

    with open("${sample_id}.plasmid_analysis.txt", "w") as sf, \\
         open("${sample_id}.plasmid.fasta", "w") as pf:

        sf.write("=" * 60 + NL)
        sf.write("PLASMID ANALYSIS REPORT" + NL)
        sf.write("=" * 60 + NL)
        sf.write(f"Sample: {sample}" + NL)
        sf.write(f"Expected CT plasmid size: ~{expected_size} bp" + NL)
        sf.write(f"Acceptable range: {min_size}-{max_size} bp" + NL + NL)

        if len(contigs) == 1 and contigs[0][0] in ["no_plasmid_detected", "assembly_failed"]:
            sf.write("STATUS: NO PLASMID DETECTED" + NL)
            pf.write(">no_plasmid" + NL + "NNNNNNNNNN" + NL)
        else:
            sf.write(f"CONTIGS FOUND: {len(contigs)}" + NL)
            sf.write("-" * 40 + NL)

            best_plasmid = None
            best_score = 0

            for header, seq, length in contigs:
                is_circular = header.split()[0] in circular_contigs or "circular" in header.lower()
                if min_size <= length <= max_size:
                    score = 100 - abs(length - expected_size) / expected_size * 100
                    if is_circular:
                        score += 20
                    if score > best_score:
                        best_score = score
                        best_plasmid = (header, seq, length, is_circular)
                circular_str = " [CIRCULAR]" if is_circular else ""
                sf.write(f"  {header.split()[0]}: {length} bp{circular_str}" + NL)

            sf.write(NL)

            if best_plasmid:
                header, seq, length, is_circular = best_plasmid
                sf.write("RESULT: CT PLASMID DETECTED" + NL)
                sf.write(f"  Length: {length} bp" + NL)
                sf.write(f"  Circular: {'Yes' if is_circular else 'No/Unknown'}" + NL)
                pf.write(f">{sample}_plasmid length={length}")
                if is_circular:
                    pf.write(" topology=circular")
                pf.write(NL)
                for i in range(0, len(seq), 80):
                    pf.write(seq[i:i+80] + NL)
            else:
                sf.write("RESULT: NO VALID CT PLASMID FOUND" + NL)
                pf.write(">no_valid_plasmid" + NL + "NNNNNNNNNN" + NL)
    """
}

/*
 * Re-orient plasmid to match reference orientation
 */
process REORIENT_PLASMID {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/${sample_id}/plasmid", mode: 'copy'

    input:
    tuple val(sample_id), path(plasmid)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid.reoriented.fasta"), emit: reoriented
    tuple val(sample_id), path("${sample_id}.orientation.txt"),          emit: orientation_log

    script:
    """
    #!/usr/bin/env python3
    import os
    import subprocess
    import re

    sample = "${sample_id}"
    plasmid_file = "${plasmid}"
    reference_file = "${params.ref_plasmid}"
    output_file = f"{sample}.plasmid.reoriented.fasta"
    log_file = f"{sample}.orientation.txt"

    def reverse_complement(seq):
        complement = dict(A='T', T='A', G='C', C='G', a='t', t='a', g='c', c='g', N='N', n='n')
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def rotate_sequence(seq, position):
        if position <= 0 or position >= len(seq):
            return seq
        return seq[position:] + seq[:position]

    def parse_fasta(filename):
        sequences = []
        with open(filename) as f:
            header = ""
            seq = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header and seq:
                        sequences.append((header, seq))
                    header = line[1:]
                    seq = ""
                else:
                    seq += line
            if header and seq:
                sequences.append((header, seq))
        return sequences

    def write_fasta(filename, header, seq):
        with open(filename, 'w') as f:
            f.write(">" + header + "\\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\\n")

    # Parse input plasmid
    plasmid_seqs = parse_fasta(plasmid_file)

    if not plasmid_seqs or "no_plasmid" in plasmid_seqs[0][0] or "no_valid" in plasmid_seqs[0][0]:
        with open(output_file, 'w') as f:
            for header, seq in (plasmid_seqs or [("no_plasmid", "NNNNNNNNNN")]):
                f.write(">" + header + "\\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\\n")
        with open(log_file, 'w') as f:
            f.write("No valid plasmid to reorient\\n")
    else:
        header, seq = plasmid_seqs[0]
        original_len = len(seq)

        # Run minimap2
        try:
            result = subprocess.run(
                ["minimap2", "-cx", "asm5", reference_file, plasmid_file],
                capture_output=True, text=True, timeout=60
            )
            paf_lines = result.stdout.strip().split("\\n")
        except Exception as e:
            write_fasta(output_file, header, seq)
            with open(log_file, 'w') as f:
                f.write("minimap2 failed: " + str(e) + "\\nKept original orientation\\n")
            paf_lines = []

        # Find best alignment by match count
        best_alignment = None
        best_matches = 0
        for line in paf_lines:
            if not line.strip():
                continue
            fields = line.split("\\t")
            if len(fields) >= 12:
                matches = int(fields[9])
                if matches > best_matches:
                    best_matches = matches
                    best_alignment = {
                        'qstart': int(fields[2]),
                        'qend': int(fields[3]),
                        'strand': fields[4],
                        'tstart': int(fields[7]),
                        'tend': int(fields[8]),
                        'matches': matches,
                        'alnlen': int(fields[10])
                    }

        if not best_alignment:
            write_fasta(output_file, header, seq)
            with open(log_file, 'w') as f:
                f.write("No alignment found\\nKept original orientation\\n")
        else:
            needs_revcomp = best_alignment['strand'] == '-'

            if needs_revcomp:
                seq = reverse_complement(seq)
                rotation_pos = original_len - best_alignment['qend']
            else:
                rotation_pos = best_alignment['qstart'] - best_alignment['tstart']

            if rotation_pos < 0:
                rotation_pos = original_len + rotation_pos
            rotation_pos = rotation_pos % original_len

            if rotation_pos > 0:
                seq = rotate_sequence(seq, rotation_pos)

            write_fasta(output_file, header + " reoriented", seq)

            with open(log_file, 'w') as f:
                f.write(f"Original length: {original_len} bp\\n")
                f.write(f"Alignment matches: {best_alignment['matches']} / {best_alignment['alnlen']}\\n")
                f.write(f"Strand: {best_alignment['strand']}\\n")
                f.write(f"Reverse complemented: {'Yes' if needs_revcomp else 'No'}\\n")
                f.write(f"Rotation applied: {rotation_pos} bp\\n")
                f.write(f"Final length: {len(seq)} bp\\n")
    """
}

/*
 * BLAST plasmid against reference plasmids
 */
process BLAST_PLASMID {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/blast:2.17.0--hb02a186_1"
    
    publishDir "${params.outdir}/${sample_id}/plasmid", mode: 'copy'

    input:
    tuple val(sample_id), path(plasmid)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid_blast.tsv"), emit: blast

    script:
    """
    if grep -q "no_plasmid\\|no_valid_plasmid\\|NNNNNNNNNN" ${plasmid}; then
        echo -e "query\\tsubject\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore" > ${sample_id}.plasmid_blast.tsv
    else
        blastn \\
            -query ${plasmid} \\
            -db ${params.plasmid_blast_db_prefix} \\
            -outfmt 6 \\
            -max_target_seqs 5 \\
            -out ${sample_id}.plasmid_blast.tsv \\
            2> blast.log
    fi
    """
}

/*
 * MLST typing of assembled plasmid
 */
process MLST_PLASMID {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pymlst:2.3.2--pyhdfd78af_0"

    publishDir "${params.outdir}/${sample_id}/plasmid", mode: 'copy'

    input:
    tuple val(sample_id), path(plasmid)

    output:
    tuple val(sample_id), path("${sample_id}.plasmid_mlst.txt"), emit: mlst

    script:
    """
    if grep -q "no_plasmid\\|no_valid_plasmid\\|NNNNNNNNNN" ${plasmid}; then
        echo -e "Sample\\tST\\tAlleles" > ${sample_id}.plasmid_mlst.txt
        echo -e "${sample_id}\\tNA\\tNA" >> ${sample_id}.plasmid_mlst.txt
    else
        claMLST search \\
            ${params.plasmid_mlst_db}/plasmid \\
            ${plasmid} > ${sample_id}.plasmid_mlst.txt 2> mlst.log || {
                echo -e "Sample\\tST\\tAlleles" > ${sample_id}.plasmid_mlst.txt
                echo -e "${sample_id}\\tNA\\tNA" >> ${sample_id}.plasmid_mlst.txt
            }
    fi
    """
}

/*
 * Collate plasmid results across all samples
 */
process COLLATE_PLASMID {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy', pattern: '*.{tsv,txt}'
    publishDir "${params.outdir}",         mode: 'copy', pattern: 'plasmid.*.all_samples.tsv'

    input:
    path ('analyses/*')
    path ('blasts/*')
    path ('mlsts/*')

    output:
    path "plasmid_summary.tsv",              emit: summary
    path "plasmid.blast.all_samples.tsv",    emit: blast_all
    path "plasmid.mlst.all_samples.tsv",     emit: mlst_all
    path "plasmid_analysis_report.txt",      emit: report

    script:
    """
    #!/usr/bin/env python3
    import csv
    import os
    from pathlib import Path

    results = []

    for analysis_file in sorted(Path("analyses").glob("*")):
        content = analysis_file.read_text()
        sample = analysis_file.stem.replace(".plasmid_analysis", "")
        row = {"sample": sample}

        if "CT PLASMID DETECTED" in content:
            row["plasmid_status"] = "DETECTED"
            for line in content.split("\\n"):
                if "Length:" in line:
                    row["plasmid_length"] = line.split(":")[1].strip().replace(" bp", "")
                if "Circular:" in line:
                    row["circular"] = line.split(":")[1].strip()
        elif "NO PLASMID DETECTED" in content:
            row["plasmid_status"] = "NOT_DETECTED"
            row["plasmid_length"] = ""
            row["circular"] = ""
        else:
            row["plasmid_status"] = "NO_VALID_PLASMID"
            row["plasmid_length"] = ""
            row["circular"] = ""

        # Parse BLAST best hit
        blast_file = Path("blasts") / f"{sample}.plasmid_blast.tsv"
        row["blast_hit"] = ""
        row["blast_identity"] = ""
        if blast_file.exists() and blast_file.stat().st_size > 0:
            with open(blast_file) as f:
                for line in f:
                    if not line.startswith("query"):
                        fields = line.strip().split("\\t")
                        if len(fields) >= 3:
                            row["blast_hit"] = fields[1]
                            row["blast_identity"] = fields[2]
                        break

        results.append(row)

    # Write summary
    with open("plasmid_summary.tsv", "w", newline="") as f:
        fieldnames = ["sample", "plasmid_status", "plasmid_length", "circular", "blast_hit", "blast_identity"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    # Collate BLAST
    with open("plasmid.blast.all_samples.tsv", "w") as f:
        f.write("query\\tsubject\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\n")
        for blast_file in sorted(Path("blasts").glob("*")):
            with open(blast_file) as bf:
                for line in bf:
                    if not line.startswith("query"):
                        f.write(line)

    # Collate MLST
    with open("plasmid.mlst.all_samples.tsv", "w") as f:
        first = True
        for mlst_file in sorted(Path("mlsts").glob("*")):
            with open(mlst_file) as mf:
                for i, line in enumerate(mf):
                    if first or i > 0:
                        f.write(line)
                first = False

    # Generate report
    detected = sum(1 for r in results if r.get("plasmid_status") == "DETECTED")
    not_detected = sum(1 for r in results if r.get("plasmid_status") == "NOT_DETECTED")

    with open("plasmid_analysis_report.txt", "w") as f:
        f.write("=" * 60 + "\\n")
        f.write("CT PLASMID ANALYSIS SUMMARY\\n")
        f.write("=" * 60 + "\\n\\n")
        f.write(f"Total samples: {len(results)}\\n")
        f.write(f"Plasmid detected: {detected}\\n")
        f.write(f"Plasmid not detected: {not_detected}\\n")
        f.write(f"Other: {len(results) - detected - not_detected}\\n\\n")
        for r in results:
            f.write(f"\\n{r['sample']}:\\n")
            f.write(f"  Status: {r.get('plasmid_status', 'UNKNOWN')}\\n")
            if r.get('plasmid_length'):
                f.write(f"  Length: {r['plasmid_length']} bp\\n")
            if r.get('blast_hit'):
                f.write(f"  Best BLAST hit: {r['blast_hit']} ({r['blast_identity']}% identity)\\n")
    """
}

/*
 * Generate per-sample plasmid PDF report
 */
process PLASMID_REPORT {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/${sample_id}/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(analysis), path(blast), path(mlst), path(mapping_stats)

    output:
    tuple val(sample_id), path("${sample_id}_plasmid_report.pdf"), emit: pdf

    script:
    """
    #!/usr/bin/env python3
    import os
    import re
    from pathlib import Path
    from datetime import datetime

    sample = "${sample_id}"

    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
        )
        HAS_REPORTLAB = True
    except ImportError:
        HAS_REPORTLAB = False

    # ---- Parse analysis file ----
    analysis_result = {"status": "UNKNOWN", "length": "", "circular": "", "deviation": ""}
    analysis_path = "${analysis}"
    if os.path.exists(analysis_path):
        content = Path(analysis_path).read_text()
        if "CT PLASMID DETECTED" in content:
            analysis_result["status"] = "DETECTED"
        elif "NO PLASMID DETECTED" in content:
            analysis_result["status"] = "NOT_DETECTED"
        elif "NO VALID CT PLASMID" in content:
            analysis_result["status"] = "NO_VALID_PLASMID"
        for line in content.split("\\n"):
            if "Length:" in line and "bp" in line:
                analysis_result["length"] = line.split(":")[1].strip()
            if "Circular:" in line:
                analysis_result["circular"] = line.split(":")[1].strip()
            if "Size deviation:" in line:
                analysis_result["deviation"] = line.split(":")[1].strip()

    # ---- Parse BLAST results ----
    blast_hits = []
    blast_path = "${blast}"
    if os.path.exists(blast_path) and os.path.getsize(blast_path) > 0:
        with open(blast_path) as f:
            for line in f:
                if line.startswith("query"):
                    continue
                fields = line.strip().split("\\t")
                if len(fields) >= 12:
                    blast_hits.append({
                        "subject": fields[1],
                        "identity": float(fields[2]),
                        "length": int(fields[3]),
                        "evalue": float(fields[10]),
                        "bitscore": float(fields[11])
                    })
        blast_hits.sort(key=lambda x: x["bitscore"], reverse=True)

    # ---- Parse mapping stats ----
    mapping = {"mapped_reads": 0, "mapped_percent": "0%"}
    mapping_path = "${mapping_stats}"
    if os.path.exists(mapping_path):
        with open(mapping_path) as f:
            for line in f:
                if "mapped (" in line and "primary" not in line:
                    parts = line.split()
                    if parts:
                        mapping["mapped_reads"] = int(parts[0])
                        if "(" in line and "%" in line:
                            pct = line.split("(")[1].split("%")[0]
                            mapping["mapped_percent"] = f"{pct}%"

    # ---- Parse MLST ----
    mlst_result = "Not available"
    mlst_path = "${mlst}"
    if os.path.exists(mlst_path) and os.path.getsize(mlst_path) > 0:
        with open(mlst_path) as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
        if len(lines) > 1:
            fields = lines[1].split("\\t")
            if len(fields) >= 2:
                mlst_result = f"ST{fields[1]}" if fields[1] != "NA" else "No match"

    # ---- Generate report ----
    if HAS_REPORTLAB:
        doc = SimpleDocTemplate(
            "${sample_id}_plasmid_report.pdf",
            pagesize=letter,
            rightMargin=0.75*inch, leftMargin=0.75*inch,
            topMargin=0.75*inch, bottomMargin=0.75*inch
        )
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle('Title', parent=styles['Heading1'],
                                     fontSize=18, alignment=TA_CENTER, spaceAfter=12)
        section_style = ParagraphStyle('Section', parent=styles['Heading2'],
                                       fontSize=14, textColor=colors.HexColor("#2C3E50"),
                                       spaceBefore=12, spaceAfter=6)
        elements = []

        elements.append(Paragraph("CT Plasmid Analysis Report", title_style))
        elements.append(Spacer(1, 0.1*inch))

        header_data = [
            ["Sample ID:", sample],
            ["Report Date:", datetime.now().strftime("%Y-%m-%d %H:%M")],
            ["Method:", "Reference-guided read extraction + Unicycler"],
        ]
        ht = Table(header_data, colWidths=[1.5*inch, 4*inch])
        ht.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ]))
        elements.append(ht)
        elements.append(Spacer(1, 0.2*inch))
        elements.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#3498DB")))
        elements.append(Spacer(1, 0.2*inch))

        # Detection result
        elements.append(Paragraph("Plasmid Detection Result", section_style))
        if analysis_result["status"] == "DETECTED":
            status_color = colors.HexColor("#27AE60")
            status_text = "PLASMID DETECTED"
        elif analysis_result["status"] == "NOT_DETECTED":
            status_color = colors.HexColor("#E74C3C")
            status_text = "NO PLASMID DETECTED"
        else:
            status_color = colors.HexColor("#F39C12")
            status_text = "UNCERTAIN"

        result_data = [
            ["Status:", status_text],
            ["Length:", analysis_result.get("length", "N/A")],
            ["Circular:", analysis_result.get("circular", "Unknown")],
        ]
        rt = Table(result_data, colWidths=[1.5*inch, 4*inch])
        rt.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 11),
            ('BACKGROUND', (1, 0), (1, 0), status_color),
            ('TEXTCOLOR', (1, 0), (1, 0), colors.white),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
        ]))
        elements.append(rt)
        elements.append(Spacer(1, 0.2*inch))

        # Mapping stats
        elements.append(Paragraph("Read Mapping to Reference Plasmid", section_style))
        mt_data = [
            ["Metric", "Value"],
            ["Reads mapped to pCT reference", f"{mapping.get('mapped_reads', 0):,}"],
            ["Mapping rate", mapping.get("mapped_percent", "0%")],
        ]
        mt = Table(mt_data, colWidths=[3*inch, 2*inch])
        mt.setStyle(TableStyle([
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#9B59B6")),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
        ]))
        elements.append(mt)
        elements.append(Spacer(1, 0.2*inch))

        # BLAST results
        elements.append(Paragraph("BLAST Results (vs 20 CT Reference Plasmids)", section_style))
        if blast_hits:
            bt_data = [["Reference Plasmid", "Identity (%)", "Alignment", "E-value"]]
            for hit in blast_hits[:5]:
                bt_data.append([
                    hit["subject"][:30],
                    f"{hit['identity']:.1f}",
                    f"{hit['length']} bp",
                    f"{hit['evalue']:.1e}"
                ])
            bt = Table(bt_data, colWidths=[2.5*inch, 1*inch, 1*inch, 1*inch])
            bt.setStyle(TableStyle([
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#27AE60")),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('BOX', (0, 0), (-1, -1), 1, colors.HexColor("#BDC3C7")),
                ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor("#BDC3C7")),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
                ('TOPPADDING', (0, 0), (-1, -1), 3),
            ]))
            elements.append(bt)
        else:
            elements.append(Paragraph("No BLAST hits found.", styles['Normal']))
        elements.append(Spacer(1, 0.2*inch))

        # MLST
        elements.append(Paragraph("Plasmid MLST", section_style))
        elements.append(Paragraph(f"<b>MLST Result:</b> {mlst_result}", styles['Normal']))
        elements.append(Spacer(1, 0.3*inch))

        # Methodology
        elements.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#BDC3C7")))
        elements.append(Spacer(1, 0.1*inch))
        method_text = (
            '<font size=8 color="#7F8C8D">'
            '<b>Methodology:</b> Reads were mapped to the reference CT plasmid (pCT, D strain, ~7,500 bp). '
            'Mapped read pairs were extracted and assembled using Unicycler. The assembly was verified '
            'by size (expected 6,000-10,000 bp) and BLAST against 20 CT reference plasmids.<br/><br/>'
            'Generated by CtGAP (Chlamydia trachomatis Genome Analysis Pipeline)'
            '</font>'
        )
        elements.append(Paragraph(method_text, styles['Normal']))
        doc.build(elements)
    else:
        # Fallback: plain text
        with open("${sample_id}_plasmid_report.pdf", "w") as f:
            f.write(f"CT Plasmid Analysis Report: {sample}\\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\\n")
            f.write("=" * 60 + "\\n\\n")
            f.write(f"Status: {analysis_result['status']}\\n")
            f.write(f"Length: {analysis_result.get('length', 'N/A')}\\n")
            f.write(f"Circular: {analysis_result.get('circular', 'Unknown')}\\n\\n")
            f.write(f"Mapped reads: {mapping.get('mapped_reads', 0):,}\\n")
            f.write(f"Mapping rate: {mapping.get('mapped_percent', '0%')}\\n\\n")
            if blast_hits:
                f.write("Top BLAST hits:\\n")
                for hit in blast_hits[:5]:
                    f.write(f"  {hit['subject']}: {hit['identity']:.1f}% identity\\n")
            f.write(f"\\nMLST: {mlst_result}\\n")
    """
}

// MAKE_PLASMID_BLASTDB removed — pre-built BLAST database used directly
// from mount point via params.plasmid_blast_db_prefix
