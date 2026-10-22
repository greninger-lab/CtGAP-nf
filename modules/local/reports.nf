/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPORTING MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules from 8-reports.smk:
      create_comprehensive_report, create_summary_csv, create_batch_summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * Per-sample comprehensive PDF report
 * Equivalent to Snakemake rule: create_comprehensive_report
 */
process COMPREHENSIVE_REPORT {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/jefffurlong/pandas-reportlab"

    publishDir "${params.outdir}/${sample_id}/reports", mode: 'copy'

    input:
    tuple val(sample_id),
          path(assembly),
          path(method_file),
          path(denovo_quast),
          path(blast_ompa),
          path(blast_secondary),
          path(mlst_chlamydiales),
          path(mlst_ct),
          path(plasmid_analysis),
          path(plasmid_blast),
          path(plasmid_mlst),
          path(fastp_json),
          path(kraken_report),
          path(coverage_tsv),
          path(annotation_txt),
          path(mixed_summary)

    output:
    tuple val(sample_id), path("${sample_id}_ctgap_report.pdf"), emit: pdf

    script:
    """
    #!/usr/bin/env python3
    import subprocess
    import sys
    import reportlab
    import json
    import csv
    import re
    import os
    from pathlib import Path
    from datetime import datetime

    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
    )
    HAS_REPORTLAB = True

    sample = "${sample_id}"
    CT_REFERENCE_SIZE = 1_043_000

    # =========================================================================
    # PARSING FUNCTIONS
    # =========================================================================
    def parse_quast(tsv_path):
        stats = {}
        if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
            return stats
        with open(tsv_path) as f:
            for line in f:
                parts = line.strip().split("\\t")
                if len(parts) >= 2:
                    key, val = parts[0], parts[1]
                    if key == "Total length":
                        stats["total_length"] = int(val.replace(",", ""))
                    elif key == "# contigs":
                        stats["num_contigs"] = int(val)
                    elif key == "N50":
                        stats["n50"] = int(val.replace(",", ""))
                    elif key == "GC (%)":
                        stats["gc_percent"] = float(val)
                    elif key == "Largest contig":
                        stats["largest_contig"] = int(val.replace(",", ""))
        return stats

    def count_gaps_in_main_contig(assembly_file):
        if not os.path.exists(assembly_file):
            return 0
        largest_seq = ""
        current_seq = ""
        with open(assembly_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if len(current_seq) > len(largest_seq):
                        largest_seq = current_seq
                    current_seq = ""
                else:
                    current_seq += line
            if len(current_seq) > len(largest_seq):
                largest_seq = current_seq
        return largest_seq.upper().count('N')

    def parse_blast(blast_file):
        results = []
        if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
            return results
        with open(blast_file) as f:
            for line in f:
                fields = line.strip().split("\\t")
                if len(fields) >= 12:
                    try:
                        identity = float(fields[2])
                    except ValueError:
                        continue
                    results.append({
                        "query": fields[0], "subject": fields[1],
                        "identity": identity, "length": int(fields[3]),
                        "evalue": float(fields[10]), "bitscore": float(fields[11])
                    })
        results.sort(key=lambda x: x["bitscore"], reverse=True)
        return results

    def parse_mlst(mlst_file):
        if not os.path.exists(mlst_file) or os.path.getsize(mlst_file) == 0:
            return "Not available"
        with open(mlst_file) as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
        if len(lines) > 1:
            fields = lines[1].split("\\t")
            if len(fields) >= 2 and fields[1] not in ["-", "NA", "ST"]:
                st = fields[1]
                return f"ST{st}" if not st.startswith("ST") else st
        return "No match found"

    def parse_fastp(json_file):
        result = {}
        if not os.path.exists(json_file) or os.path.getsize(json_file) == 0:
            return result
        try:
            with open(json_file) as f:
                data = json.load(f)
            if "summary" in data:
                before = data["summary"].get("before_filtering", {})
                after = data["summary"].get("after_filtering", {})
                result["raw_reads"] = before.get("total_reads", 0)
                result["filtered_reads"] = after.get("total_reads", 0)
                result["q30_rate"] = after.get("q30_rate", 0)
                result["gc_content"] = after.get("gc_content", 0)
        except Exception:
            pass
        return result

    def parse_kraken(report_file):
        results = []
        if not os.path.exists(report_file) or os.path.getsize(report_file) == 0:
            return results
        with open(report_file) as f:
            for line in f:
                parts = line.strip().split("\\t")
                if len(parts) >= 6:
                    try:
                        pct = float(parts[0])
                    except ValueError:
                        continue
                    name = parts[5].strip()
                    if pct >= 1.0 and name:
                        results.append({"name": name, "percent": pct})
        return sorted(results, key=lambda x: x["percent"], reverse=True)[:10]

    def parse_coverage(cov_file):
        result = {}
        if not os.path.exists(cov_file) or os.path.getsize(cov_file) == 0:
            return result
        with open(cov_file) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\\t")
                if len(parts) >= 8:
                    result["mean_depth"] = float(parts[7]) if parts[7] != "nan" else 0
                    result["coverage_pct"] = float(parts[6]) if parts[6] != "nan" else 0
                    break
        return result

    def parse_plasmid_analysis(analysis_file):
        result = {"status": "Unknown", "length": "", "circular": ""}
        if not os.path.exists(analysis_file):
            return result
        content = Path(analysis_file).read_text()
        if "CT PLASMID DETECTED" in content:
            result["status"] = "DETECTED"
        elif "NO PLASMID DETECTED" in content:
            result["status"] = "NOT_DETECTED"
        elif "NO VALID CT PLASMID" in content:
            result["status"] = "NO_VALID"
        for line in content.split("\\n"):
            if "Length:" in line and "bp" in line:
                m = re.search(r"Length:\\s*(\\d+)", line)
                if m:
                    result["length"] = m.group(1)
            if "Circular:" in line:
                result["circular"] = line.split(":")[1].strip()
        return result

    def parse_annotation(annotation_txt):
        if not os.path.exists(annotation_txt):
            return None
        content = Path(annotation_txt).read_text()
        result = {}
        patterns = {
            "cds": r"CDSs:\\s*(\\d+)", "trnas": r"tRNAs:\\s*(\\d+)",
            "rrnas": r"rRNAs:\\s*(\\d+)", "ncrnas": r"ncRNAs:\\s*(\\d+)",
            "hypotheticals": r"hypotheticals:\\s*(\\d+)", "pseudogenes": r"pseudogenes:\\s*(\\d+)",
        }
        for key, pattern in patterns.items():
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                result[key] = match.group(1)
        return result if result else None

    def parse_mixed_detection(summary_file):
        if not os.path.exists(summary_file):
            return None
        content = Path(summary_file).read_text()
        result = {}
        m = re.search(r"Heterozygous-like sites:\\s*(\\d+)", content)
        if m:
            result["het_sites"] = int(m.group(1))
        m = re.search(r"Mean allele frequency:\\s*([\\d.]+|N/A)", content)
        if m:
            result["mean_af"] = m.group(1)
        m = re.search(r"Main contig analyzed:\\s*(\\S+)\\s*\\((\\d+)\\s*bp\\)", content)
        if m:
            result["main_contig_size"] = int(m.group(2))
        m = re.search(r"Position range:\\s*(\\d+)", content)
        if m:
            result["position_range"] = int(m.group(1))
        if "GENOME-WIDE" in content:
            result["distribution"] = "GENOME-WIDE"
        elif "LOCALIZED" in content:
            result["distribution"] = "LOCALIZED"
        if "STATUS: MIXED_STRAIN_WARNING" in content:
            result["status"] = "MIXED_WARNING"
        elif "STATUS: LOW_HETEROZYGOSITY" in content:
            result["status"] = "LOW_HET"
        elif "STATUS: SINGLE_STRAIN" in content:
            result["status"] = "SINGLE"
        return result if result else None

    # =========================================================================
    # COLLECT DATA
    # =========================================================================
    method = Path("${method_file}").read_text().strip() if os.path.exists("${method_file}") else "unknown"
    quast = parse_quast("${denovo_quast}")
    quast["main_contig_gaps"] = count_gaps_in_main_contig("${assembly}")
    blast_ompa = parse_blast("${blast_ompa}")
    blast_sec = parse_blast("${blast_secondary}")
    mlst_chlam = parse_mlst("${mlst_chlamydiales}")
    mlst_ct = parse_mlst("${mlst_ct}")
    fastp = parse_fastp("${fastp_json}")
    kraken = parse_kraken("${kraken_report}")
    coverage = parse_coverage("${coverage_tsv}")
    plasmid = parse_plasmid_analysis("${plasmid_analysis}")
    plasmid_blast = parse_blast("${plasmid_blast}")
    plasmid_mlst = parse_mlst("${plasmid_mlst}")

    # Optional sections - read from staged input files
    annotation = parse_annotation("${annotation_txt}")
    mixed = parse_mixed_detection("${mixed_summary}")

    # Extract genovar
    genovar = "Unknown"
    genovar_identity = "N/A"
    if blast_ompa:
        subject = blast_ompa[0].get("subject", "")
        m = re.match(r"([A-Za-z]+)_", subject)
        if m:
            genovar = m.group(1)
        genovar_identity = f"{blast_ompa[0].get('identity', 0):.1f}%"

    # Plasmid summary
    plasmid_summary = plasmid.get("status", "Unknown")
    if plasmid_summary == "DETECTED" and plasmid.get("length"):
        plasmid_summary = f"DETECTED ({plasmid['length']} bp)"

    # =========================================================================
    # GENERATE REPORT
    # =========================================================================
    if HAS_REPORTLAB:
        doc = SimpleDocTemplate("${sample_id}_ctgap_report.pdf", pagesize=letter,
                                rightMargin=0.75*inch, leftMargin=0.75*inch,
                                topMargin=0.75*inch, bottomMargin=0.75*inch)
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle('Title', parent=styles['Heading1'], fontSize=18, alignment=TA_CENTER, spaceAfter=12)
        section_style = ParagraphStyle('Section', parent=styles['Heading2'], fontSize=14,
                                       textColor=colors.HexColor("#2C3E50"), spaceBefore=12, spaceAfter=6)
        subsection_style = ParagraphStyle('Subsection', parent=styles['Heading3'], fontSize=11,
                                          textColor=colors.HexColor("#34495E"), spaceBefore=8, spaceAfter=4)
        kv_style = ParagraphStyle('KV', parent=styles['Normal'], fontSize=10, leftIndent=20)

        elements = []

        # Title
        elements.append(Paragraph("CtGAP-nf Analysis Report", title_style))
        elements.append(Spacer(1, 0.1*inch))
        ht = Table([["Sample ID:", sample], ["Report Date:", datetime.now().strftime("%Y-%m-%d %H:%M")],
                     ["Pipeline:", "CtGAP-nf v1.0"]], colWidths=[1.5*inch, 4*inch])
        ht.setStyle(TableStyle([('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'), ('FONTSIZE', (0,0), (-1,-1), 10),
                                 ('BOTTOMPADDING', (0,0), (-1,-1), 4)]))
        elements.append(ht)
        elements.append(Spacer(1, 0.2*inch))
        elements.append(HRFlowable(width="100%", thickness=2, color=colors.HexColor("#3498DB")))
        elements.append(Spacer(1, 0.2*inch))

        # Key Findings
        elements.append(Paragraph("Key Findings", section_style))
        fd = [["ompA Genotype:", genovar, "Identity:", genovar_identity],
              ["Assembly Method:", method, "Contigs:", str(quast.get("num_contigs", "N/A"))],
              ["Plasmid:", plasmid_summary, "N50:", f"{quast.get('n50', 0):,} bp"]]
        ft = Table(fd, colWidths=[1.5*inch, 1.5*inch, 1.3*inch, 1.7*inch])
        ft.setStyle(TableStyle([('FONTNAME', (0,0), (0,-1), 'Helvetica-Bold'), ('FONTNAME', (2,0), (2,-1), 'Helvetica-Bold'),
                                 ('FONTSIZE', (0,0), (-1,-1), 10), ('BACKGROUND', (0,0), (-1,-1), colors.HexColor("#ECF0F1")),
                                 ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                 ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7")),
                                 ('BOTTOMPADDING', (0,0), (-1,-1), 6), ('TOPPADDING', (0,0), (-1,-1), 6)]))
        elements.append(ft)
        elements.append(Spacer(1, 0.3*inch))

        # 1. Assembly
        elements.append(Paragraph("1. Assembly Analysis", section_style))
        elements.append(Paragraph(f"<b>Selected Method:</b> {method}", kv_style))
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph("Assembly Quality Metrics (QUAST)", subsection_style))
        if quast:
            gap_count = quast.get("main_contig_gaps", 0)
            qd = [["Metric", "Value"], ["Total Length", f"{quast.get('total_length', 0):,} bp"],
                  ["Number of Contigs", str(quast.get("num_contigs", "N/A"))],
                  ["N50", f"{quast.get('n50', 0):,} bp"], ["Gaps in Main Contig", f"{gap_count:,} bp"],
                  ["GC Content", f"{quast.get('gc_percent', 0):.2f}%"]]
            qt = Table(qd, colWidths=[2.5*inch, 2.5*inch])
            qt.setStyle(TableStyle([('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                     ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#3498DB")),
                                     ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('FONTSIZE', (0,0), (-1,-1), 9),
                                     ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                     ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7")),
                                     ('ALIGN', (1,1), (1,-1), 'RIGHT'), ('BOTTOMPADDING', (0,0), (-1,-1), 4)]))
            elements.append(qt)
        elements.append(Spacer(1, 0.2*inch))

        # 2. Typing
        elements.append(Paragraph("2. Typing Results", section_style))
        elements.append(Paragraph("ompA Genotype (BLAST)", subsection_style))
        if blast_ompa:
            top = blast_ompa[0]
            elements.append(Paragraph(f"<b>Best Match:</b> {top['subject']}<br/><b>Identity:</b> {top['identity']:.2f}%<br/>"
                                      f"<b>Alignment Length:</b> {top['length']} bp<br/><b>Genotype:</b> {genovar}", kv_style))
        else:
            elements.append(Paragraph("No BLAST hits found", kv_style))
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph("MLST", subsection_style))
        elements.append(Paragraph(f"<b>Chlamydiales:</b> {mlst_chlam}<br/><b>C. trachomatis:</b> {mlst_ct}", kv_style))
        elements.append(Spacer(1, 0.2*inch))

        # 3. Plasmid
        elements.append(Paragraph("3. Plasmid Analysis", section_style))
        status_color = "#27AE60" if plasmid["status"] == "DETECTED" else "#E74C3C" if plasmid["status"] == "NOT_DETECTED" else "#F39C12"
        status_text = "PLASMID DETECTED" if plasmid["status"] == "DETECTED" else "NO PLASMID DETECTED" if plasmid["status"] == "NOT_DETECTED" else plasmid["status"]
        elements.append(Paragraph(f'<b>Status:</b> <font color="{status_color}"><b>{status_text}</b></font><br/>'
                                  f'<b>Length:</b> {plasmid.get("length", "N/A")} bp<br/>'
                                  f'<b>Circular:</b> {plasmid.get("circular", "Unknown")}', kv_style))
        if plasmid_blast:
            elements.append(Spacer(1, 0.1*inch))
            elements.append(Paragraph("Plasmid BLAST Results", subsection_style))
            pbd = [["Subject", "Identity", "Length", "E-value"]]
            for hit in plasmid_blast[:5]:
                pbd.append([hit["subject"][:30], f"{hit['identity']:.1f}%", str(hit["length"]), f"{hit['evalue']:.1e}"])
            pbt = Table(pbd, colWidths=[2.5*inch, 1*inch, 1*inch, 1.2*inch])
            pbt.setStyle(TableStyle([('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                      ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#27AE60")),
                                      ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('FONTSIZE', (0,0), (-1,-1), 8),
                                      ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                      ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7"))]))
            elements.append(pbt)
        if plasmid_mlst and plasmid_mlst != "Not available":
            elements.append(Paragraph(f"<b>Plasmid MLST:</b> {plasmid_mlst}", kv_style))
        elements.append(Spacer(1, 0.2*inch))

        # 4. Annotation (optional)
        sec_num = 4
        if annotation:
            elements.append(Paragraph(f"{sec_num}. Genome Annotation (Bakta)", section_style))
            cds = int(annotation.get("cds", 0))
            hypo = int(annotation.get("hypotheticals", 0))
            annotated = cds - hypo if cds > hypo else cds
            rate = (annotated / cds * 100) if cds > 0 else 0
            ad = [["Feature", "Count"], ["CDS", annotation.get("cds", "N/A")],
                  ["tRNAs", annotation.get("trnas", "N/A")], ["rRNAs", annotation.get("rrnas", "N/A")],
                  ["ncRNAs", annotation.get("ncrnas", "N/A")], ["Hypotheticals", annotation.get("hypotheticals", "N/A")],
                  ["Pseudogenes", annotation.get("pseudogenes", "N/A")]]
            at = Table(ad, colWidths=[3*inch, 2*inch])
            at.setStyle(TableStyle([('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                     ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#16A085")),
                                     ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('FONTSIZE', (0,0), (-1,-1), 9),
                                     ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                     ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7")),
                                     ('ALIGN', (1,1), (1,-1), 'RIGHT')]))
            elements.append(at)
            elements.append(Paragraph(f"<b>Annotation Efficiency:</b> {rate:.1f}% of CDS have functional annotation", kv_style))
            elements.append(Spacer(1, 0.2*inch))
            sec_num += 1

        # QC Section
        elements.append(Paragraph(f"{sec_num}. Quality Control", section_style))
        if fastp:
            elements.append(Paragraph("Read Statistics", subsection_style))
            raw = fastp.get("raw_reads", 0)
            filt = fastp.get("filtered_reads", 0)
            retention = 100 * filt / max(raw, 1)
            qcd = [["Metric", "Value"], ["Raw Reads", f"{raw:,}"], ["Filtered Reads", f"{filt:,}"],
                   ["Retention Rate", f"{retention:.1f}%"], ["Q30 Rate", f"{100*fastp.get('q30_rate', 0):.1f}%"],
                   ["GC Content", f"{100*fastp.get('gc_content', 0):.1f}%"]]
            qct = Table(qcd, colWidths=[2.5*inch, 2.5*inch])
            qct.setStyle(TableStyle([('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                      ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#9B59B6")),
                                      ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('FONTSIZE', (0,0), (-1,-1), 9),
                                      ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                      ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7")),
                                      ('ALIGN', (1,1), (1,-1), 'RIGHT')]))
            elements.append(qct)

        if coverage or quast:
            elements.append(Spacer(1, 0.1*inch))
            elements.append(Paragraph("Coverage Statistics", subsection_style))
            n50 = quast.get("n50", 0)
            genome_cov = (n50 / CT_REFERENCE_SIZE * 100) if n50 > 0 else 0
            elements.append(Paragraph(f"<b>Mean Depth:</b> {coverage.get('mean_depth', 0):.1f}x<br/>"
                                      f"<b>Contiguous Assembly (N50):</b> {n50:,} bp<br/>"
                                      f"<b>Genome Coverage:</b> {genome_cov:.1f}%*", kv_style))
            if n50 > 1_060_000:
                elements.append(Paragraph(f'<font color="#E74C3C"><b>⚠ QC FLAG:</b> N50 ({n50:,} bp) exceeds expected Ct genome size (&gt;1.06 Mb).</font>', kv_style))
            elements.append(Paragraph('<font size=8 color="#7F8C8D">*N50 relative to standard C. trachomatis genome size (1.043 Mb)</font>', styles['Normal']))

        if kraken:
            elements.append(Spacer(1, 0.1*inch))
            elements.append(Paragraph("Kraken2 Classification (Top Hits)", subsection_style))
            kd = [["Taxon", "Percentage"]]
            for hit in kraken[:5]:
                kd.append([hit["name"][:40], f"{hit['percent']:.1f}%"])
            kt = Table(kd, colWidths=[4*inch, 1.5*inch])
            kt.setStyle(TableStyle([('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                     ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#E67E22")),
                                     ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('FONTSIZE', (0,0), (-1,-1), 8),
                                     ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                     ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7")),
                                     ('ALIGN', (1,1), (1,-1), 'RIGHT')]))
            elements.append(kt)
        sec_num += 1

        # Mixed Strain Detection (optional)
        if mixed:
            elements.append(Spacer(1, 0.2*inch))
            elements.append(Paragraph(f"{sec_num}. Mixed Strain Analysis", section_style))
            het = mixed.get("het_sites", 0)
            mean_af = mixed.get("mean_af", "N/A")
            status = mixed.get("status", "Unknown")
            pos_range = mixed.get("position_range", 0)
            dist = mixed.get("distribution", "N/A")
            range_str = f"{pos_range/1000:.1f} kb" if pos_range > 1000 else f"{pos_range} bp"

            md = [["Metric", "Value"],
                  ["Heterozygous Sites", str(het)],
                  ["Mean Allele Frequency", str(mean_af)],
                  ["Position Spread", range_str]]
            mt = Table(md, colWidths=[2.5*inch, 2.5*inch])
            mt.setStyle(TableStyle([('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
                                     ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#8E44AD")),
                                     ('TEXTCOLOR', (0,0), (-1,0), colors.white), ('FONTSIZE', (0,0), (-1,-1), 9),
                                     ('BOX', (0,0), (-1,-1), 1, colors.HexColor("#BDC3C7")),
                                     ('INNERGRID', (0,0), (-1,-1), 0.5, colors.HexColor("#BDC3C7"))]))
            elements.append(mt)
            elements.append(Spacer(1, 0.1*inch))

            if status == "MIXED_WARNING":
                dist_text = ""
                if dist == "GENOME-WIDE":
                    dist_text = "<br/><br/><b>Distribution:</b> <font color='#27AE60'>GENOME-WIDE</font> - Variants are spread across the genome, consistent with true mixed infection."
                elif dist == "LOCALIZED":
                    dist_text = "<br/><br/><b>Distribution:</b> <font color='#F39C12'>LOCALIZED</font> - Variants are clustered in one region, could indicate recombination or localized artifact."

                elements.append(Paragraph(
                    f'<b>Status:</b> <font color="#E74C3C">POSSIBLE MIXED INFECTION</font><br/><br/>'
                    f'<b>Interpretation:</b> This sample shows {het} genomic positions where sequencing reads '
                    f'support two different alleles at intermediate frequencies. Since <i>C. trachomatis</i> is haploid, '
                    f'this pattern suggests the presence of multiple distinct strains in the sample.{dist_text}<br/><br/>'
                    f'<b>Recommendations:</b><br/>'
                    f'&bull; Interpret typing results (genotype, MLST, Ct strain) with caution - they may represent a consensus<br/>'
                    f'&bull; Consider that plasmid detection may be affected by strain mixture<br/>'
                    f'&bull; Review clinical context for possible co-infection or reinfection', kv_style))
            elif status in ("SINGLE", "LOW_HET"):
                elements.append(Paragraph(
                    f'<b>Status:</b> <font color="#27AE60">SINGLE STRAIN (LOW HETEROZYGOSITY)</font><br/><br/>'
                    f'<b>Interpretation:</b> This sample shows minimal heterozygosity ({het} sites), '
                    f'consistent with a single <i>C. trachomatis</i> strain. Typing results can be interpreted '
                    f'with normal confidence.', kv_style))
            else:
                elements.append(Paragraph(
                    f'<b>Status:</b> <font color="#F39C12">INCONCLUSIVE</font><br/><br/>'
                    f'<b>Interpretation:</b> Mixed strain analysis could not determine a clear result.', kv_style))

            # Methodology note
            elements.append(Spacer(1, 0.1*inch))
            method_note_style = ParagraphStyle('MethodNote', parent=styles['Normal'], fontSize=8,
                                                textColor=colors.HexColor("#7F8C8D"), leftIndent=20)
            elements.append(Paragraph(
                '<b>Method:</b> Reads were mapped to the sample\'s main chromosome (largest contig) and variants called. '
                'Heterozygous-like sites are positions where 15-85% of reads support an alternate allele. '
                '<b>Quality filters:</b> Base quality &ge;20, minimum 2 reads on each strand (reduces sequencing errors '
                'and strand bias artifacts). Analysis focuses on main chromosome only (excludes plasmid and small fragments). '
                'In a haploid organism, sites passing these filters indicate multiple genomes.<br/><br/>'
                'Generated by CtGAP-nf', method_note_style))

        # Footer
        elements.append(Spacer(1, 0.3*inch))
        elements.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor("#BDC3C7")))
        elements.append(Spacer(1, 0.1*inch))
        elements.append(Paragraph(f'<font size=8 color="#7F8C8D">Generated by CtGAP-nf | {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</font>',
                                  ParagraphStyle('Footer', alignment=TA_CENTER)))
        doc.build(elements)

    else:
        # Fallback (should not reach here - reportlab is installed above)
        with open("${sample_id}_ctgap_report.pdf", "w") as f:
            f.write(f"CtGAP-nf Report: {sample} - reportlab installation failed\\n")
    """
}

/*
 * Master CSV summarizing all samples
 * Equivalent to Snakemake rule: create_summary_csv
 */
process SUMMARY_CSV {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('methods/*')
    path ('quast_reports/*')
    path ('blast_tabs/*')
    path ('mlst_ct/*')
    path ('plasmid_analyses/*')

    output:
    path "ctgap_all_samples_summary.csv", emit: csv

    script:
    """
    #!/usr/bin/env python3
    import csv
    import re
    import os
    from pathlib import Path

    rows = []

    # Index all files by sample name for efficient lookup
    def find_file_for_sample(directory, sample, extensions=None):
        # Find a file matching sample name in directory
        for f in Path(directory).glob("*"):
            fname = f.name
            if sample in fname:
                if extensions is None or any(fname.endswith(ext) for ext in extensions):
                    return f
        return None

    # Build sample list from method files
    for method_file in sorted(Path("methods").glob("*.method.txt")):
        sample = method_file.stem.replace(".method", "")
        row = {"sample": sample}

        # 1. Assembly method
        row["assembly_method"] = method_file.read_text().strip()

        # 2. QUAST metrics - find report matching this sample
        row["total_length"] = ""
        row["num_contigs"] = ""
        row["n50"] = ""
        row["gc_percent"] = ""
        qf = find_file_for_sample("quast_reports", sample, [".tsv"])
        if qf is None:
            # Fallback: if only one report per sample, files may just be report.tsv
            # Try all files
            for candidate in Path("quast_reports").glob("*"):
                qf = candidate
                break
        if qf and qf.exists():
            with open(qf) as f:
                for line in f:
                    parts = line.strip().split("\\t")
                    if len(parts) >= 2:
                        key, val = parts[0], parts[1]
                        if key == "Total length":
                            row["total_length"] = val.replace(",", "")
                        elif key == "# contigs":
                            row["num_contigs"] = val
                        elif key == "N50":
                            row["n50"] = val.replace(",", "")
                        elif key == "GC (%)":
                            row["gc_percent"] = val

        # 3. ompA genotype - find BLAST tab matching this sample
        row["genotype"] = ""
        row["genotype_identity"] = ""
        bf = find_file_for_sample("blast_tabs", sample, [".tab"])
        if bf and bf.exists() and bf.stat().st_size > 0:
            with open(bf) as f:
                lines = [l for l in f.readlines() if l.strip()]
                if lines:
                    best = sorted(lines, key=lambda x: float(x.split("\\t")[11]) if len(x.split("\\t")) > 11 else 0, reverse=True)[0]
                    fields = best.strip().split("\\t")
                    if len(fields) >= 3:
                        m = re.match(r"([A-Za-z]+)_", fields[1])
                        if m:
                            row["genotype"] = m.group(1)
                        row["genotype_identity"] = fields[2]

        # 4. MLST - find MLST result matching this sample
        row["mlst_st"] = ""
        mf = find_file_for_sample("mlst_ct", sample, [".txt"])
        if mf and mf.exists() and mf.stat().st_size > 0:
            with open(mf) as f:
                lines = [l.strip() for l in f if l.strip()]
                if len(lines) > 1:
                    fields = lines[1].split("\\t")
                    if len(fields) >= 2 and fields[1] not in ["-", "NA", "ST"]:
                        row["mlst_st"] = f"ST{fields[1]}" if not fields[1].startswith("ST") else fields[1]

        # 5. Plasmid - find analysis matching this sample
        row["plasmid_status"] = "Unknown"
        row["plasmid_length"] = ""
        pf = find_file_for_sample("plasmid_analyses", sample, [".txt"])
        if pf and pf.exists():
            content = pf.read_text()
            if "CT PLASMID DETECTED" in content:
                row["plasmid_status"] = "DETECTED"
                m = re.search(r"Length:\\s*(\\d+)", content)
                if m:
                    row["plasmid_length"] = m.group(1)
            elif "NO PLASMID DETECTED" in content:
                row["plasmid_status"] = "NOT_DETECTED"
            elif "NO VALID CT PLASMID" in content:
                row["plasmid_status"] = "NO_VALID"

        rows.append(row)

    fieldnames = [
        "sample", "genotype", "genotype_identity", "mlst_st",
        "assembly_method", "total_length", "num_contigs", "n50", "gc_percent",
        "plasmid_status", "plasmid_length"
    ]

    with open("ctgap_all_samples_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda x: x["sample"]))
    """
}

/*
 * Batch summary text report
 * Equivalent to Snakemake rule: create_batch_summary
 */
process BATCH_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path summary_csv

    output:
    path "ctgap_batch_summary.txt", emit: report

    script:
    """
    total=\$(tail -n +2 ${summary_csv} | wc -l | tr -d ' ')

    cat > ctgap_batch_summary.txt << ENDHEADER
======================================================================
CtGAP-nf BATCH ANALYSIS SUMMARY
======================================================================
Generated: \$(date '+%Y-%m-%d %H:%M:%S')
Total samples: \$total

----------------------------------------
GENOTYPE DISTRIBUTION
----------------------------------------
ENDHEADER

    # Genotype distribution (column 2)
    tail -n +2 ${summary_csv} | awk -F',' '{print \$2}' | grep -v '^\$' | sort | uniq -c | sort -rn | while read count geno; do
        pct=\$(awk "BEGIN {printf \\"%.1f\\", (\$count / \$total) * 100}")
        echo "  \$geno: \$count samples (\${pct}%)" >> ctgap_batch_summary.txt
    done

    cat >> ctgap_batch_summary.txt << ENDMETHOD

----------------------------------------
ASSEMBLY METHOD SELECTION
----------------------------------------
ENDMETHOD

    # Assembly method distribution (column 5)
    tail -n +2 ${summary_csv} | awk -F',' '{print \$5}' | grep -v '^\$' | sort | uniq -c | sort -rn | while read count method; do
        pct=\$(awk "BEGIN {printf \\"%.1f\\", (\$count / \$total) * 100}")
        echo "  \$method: \$count samples (\${pct}%)" >> ctgap_batch_summary.txt
    done

    cat >> ctgap_batch_summary.txt << ENDPLASMID

----------------------------------------
PLASMID DETECTION
----------------------------------------
ENDPLASMID

    # Plasmid status distribution (column 10)
    tail -n +2 ${summary_csv} | awk -F',' '{print \$10}' | grep -v '^\$' | sort | uniq -c | sort -rn | while read count status; do
        pct=\$(awk "BEGIN {printf \\"%.1f\\", (\$count / \$total) * 100}")
        echo "  \$status: \$count samples (\${pct}%)" >> ctgap_batch_summary.txt
    done

    cat >> ctgap_batch_summary.txt << ENDQUALITY

----------------------------------------
ASSEMBLY QUALITY METRICS
----------------------------------------
ENDQUALITY

    # Total length stats (column 6)
    lengths=\$(tail -n +2 ${summary_csv} | awk -F',' '\$6 ~ /^[0-9]+\$/ {print \$6}' | sort -n)
    if [ -n "\$lengths" ]; then
        len_min=\$(echo "\$lengths" | head -1)
        len_max=\$(echo "\$lengths" | tail -1)
        len_count=\$(echo "\$lengths" | wc -l | tr -d ' ')
        len_median=\$(echo "\$lengths" | sed -n "\$((len_count / 2 + 1))p")
        echo "  Total length: \$len_min - \$len_max bp (median: \$len_median)" >> ctgap_batch_summary.txt
    fi

    # N50 stats (column 8)
    n50s=\$(tail -n +2 ${summary_csv} | awk -F',' '\$8 ~ /^[0-9]+\$/ {print \$8}' | sort -n)
    if [ -n "\$n50s" ]; then
        n50_min=\$(echo "\$n50s" | head -1)
        n50_max=\$(echo "\$n50s" | tail -1)
        n50_count=\$(echo "\$n50s" | wc -l | tr -d ' ')
        n50_median=\$(echo "\$n50s" | sed -n "\$((n50_count / 2 + 1))p")
        echo "  N50: \$n50_min - \$n50_max bp (median: \$n50_median)" >> ctgap_batch_summary.txt
    fi

    # Contig count stats (column 7)
    contigs=\$(tail -n +2 ${summary_csv} | awk -F',' '\$7 ~ /^[0-9]+\$/ {print \$7}' | sort -n)
    if [ -n "\$contigs" ]; then
        c_min=\$(echo "\$contigs" | head -1)
        c_max=\$(echo "\$contigs" | tail -1)
        c_count=\$(echo "\$contigs" | wc -l | tr -d ' ')
        c_median=\$(echo "\$contigs" | sed -n "\$((c_count / 2 + 1))p")
        echo "  Contigs: \$c_min - \$c_max (median: \$c_median)" >> ctgap_batch_summary.txt
    fi

    cat >> ctgap_batch_summary.txt << ENDFOOTER

======================================================================
For detailed per-sample reports, see individual PDF reports.
For full data, see: ctgap_all_samples_summary.csv
======================================================================
ENDFOOTER
    """
}
