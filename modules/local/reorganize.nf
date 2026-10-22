/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REORGANIZE - Create clean user-facing output structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: reorganize (9-reorganize.smk)
    Creates 01-04 directory structure with symlinks/copies to final outputs.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process REORGANIZE {
    tag "final_output"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    // Summary reports
    path summary_csv
    path batch_summary
    // Per-sample PDFs
    path ('pdfs/*')
    // Per-sample assemblies
    path ('assemblies/*')
    // Phylogeny
    path tree
    // Plasmid
    path ('plasmid_fastas/*')
    path ('plasmid_pdfs/*')
    path plasmid_summary
    path plasmid_blast_all
    path plasmid_mlst_all
    path plasmid_report
    // Typing collated results
    path ompa_blast_tsv
    path secondary_blast_tsv
    path mlst_generic_tsv
    path mlst_ct_tsv
    // Coverage
    path ('coverage_files/*')
    // Assembly selection
    path assembly_methods_summary
    path ('method_files/*')
    // MultiQC
    path multiqc_report
    // Scaffold summaries
    path ('scaffold_summaries/*')
    // Alignment summaries
    path ('alignment_summaries/*')
    // Annotation (optional)
    path annotation_summary
    path annotation_report
    path annotation_coverage_check
    // Mixed detection (optional)
    path mixed_summary
    path mixed_report

    output:
    path "00_README.txt",       emit: readme
    path "01_per_sample/**",    emit: per_sample
    path "02_reports/**",       emit: reports
    path "03_results/**",       emit: results
    path "04_phylogeny/**",     emit: phylogeny
    path "05_annotation/**",    emit: annotation, optional: true

    script:
    def org = params.organism
    """
#!/usr/bin/env python3
import os
import shutil
from pathlib import Path

org = "${org}"

# =========================================================
# Create directory structure matching 9-reorganize.smk
# =========================================================
dirs = {
    "per_sample":  Path("01_per_sample"),
    "reports":     Path("02_reports"),
    "results":     Path("03_results"),
    "phylogeny":   Path("04_phylogeny"),
}

for d in dirs.values():
    d.mkdir(parents=True, exist_ok=True)

# Subdirectories
(dirs["reports"] / "per_sample").mkdir(exist_ok=True)
(dirs["reports"] / "summaries").mkdir(exist_ok=True)
(dirs["reports"] / "qc").mkdir(exist_ok=True)
(dirs["results"] / "typing").mkdir(exist_ok=True)
(dirs["results"] / "coverage").mkdir(exist_ok=True)
(dirs["results"] / "plasmid").mkdir(exist_ok=True)
(dirs["results"] / "plasmid" / "per_sample").mkdir(exist_ok=True)

def safe_copy(src, dst):
    src = Path(src)
    dst = Path(dst)
    if src.exists() and src.stat().st_size > 0:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(src), str(dst))
        return True
    return False

# =========================================================
# 01_per_sample: Per-sample assemblies and reports
# =========================================================
# Assemblies - extract sample name from the assembly filename
# Nextflow stages files with their original names from publishDir
# Possible patterns: assembly.filtered.fasta, assembly.fasta,
# denovo_{sample}.final.fasta, ref-denovo_{sample}.final.fasta
for asm in sorted(Path("assemblies").glob("*")):
    name = asm.stem  # e.g. "assembly.filtered" or "denovo_SRR30151523.final"
    
    # Try to extract sample name from known patterns
    sample_name = None
    
    # Pattern: denovo_{sample}.final or ref-denovo_{sample}.final
    import re
    m = re.match(r"(?:denovo|ref-denovo)_(.+?)\.final", name)
    if m:
        sample_name = m.group(1)
    
    # Pattern: assembly.filtered or assembly (no sample name embedded)
    # In this case, we need to get sample name from other files (PDFs, method files)
    if not sample_name and name.startswith("assembly"):
        # Look at PDF files to find sample names
        pdf_samples = set()
        for pdf in Path("pdfs").glob("*_ctgap_report.pdf"):
            pdf_samples.add(pdf.stem.replace("_ctgap_report", ""))
        if len(pdf_samples) == 1:
            sample_name = pdf_samples.pop()
        elif len(pdf_samples) > 1:
            # Multiple samples - try method files
            for mf in Path("method_files").glob("*.method.txt"):
                sample_name = mf.stem.replace(".method", "")
                break
    
    if not sample_name:
        # Last resort: use the filename as-is but clean it up
        sample_name = name.replace("assembly.filtered", "unknown").replace("assembly", "unknown")
    
    sample_dir = dirs["per_sample"] / sample_name
    sample_dir.mkdir(exist_ok=True)
    safe_copy(asm, sample_dir / "assembly.fasta")

# PDF reports
for pdf in sorted(Path("pdfs").glob("*_ctgap_report.pdf")):
    sample = pdf.stem.replace("_ctgap_report", "")
    sample_dir = dirs["per_sample"] / sample
    sample_dir.mkdir(exist_ok=True)
    safe_copy(pdf, sample_dir / pdf.name)
    # Also copy to reports/per_sample
    safe_copy(pdf, dirs["reports"] / "per_sample" / pdf.name)

# Plasmid PDFs per sample
for pdf in sorted(Path("plasmid_pdfs").glob("*_plasmid_report.pdf")):
    sample = pdf.stem.replace("_plasmid_report", "")
    sample_dir = dirs["per_sample"] / sample
    if sample_dir.exists():
        safe_copy(pdf, sample_dir / pdf.name)
    # Also to reports/per_sample
    safe_copy(pdf, dirs["reports"] / "per_sample" / pdf.name)

# Plasmid fastas per sample
for pf in sorted(Path("plasmid_fastas").glob("*")):
    # Extract sample from filename like {sample}.plasmid.fasta
    name = pf.name
    sample = name.split(".plasmid")[0] if ".plasmid" in name else pf.stem
    sample_dir = dirs["per_sample"] / sample
    if sample_dir.exists():
        safe_copy(pf, sample_dir / pf.name)
    # Also to results/plasmid/per_sample
    safe_copy(pf, dirs["results"] / "plasmid" / "per_sample" / pf.name)

# =========================================================
# 02_reports: Summary reports
# =========================================================
# Main summaries
safe_copy("${summary_csv}", dirs["reports"] / "summaries" / "ctgap_all_samples.csv")
safe_copy("${batch_summary}", dirs["reports"] / "summaries" / "batch_summary.txt")
safe_copy("${assembly_methods_summary}", dirs["reports"] / "summaries" / "assembly_methods.txt")

# Scaffold summaries
for sf in sorted(Path("scaffold_summaries").glob("*")):
    safe_copy(sf, dirs["reports"] / "summaries" / sf.name)

# Alignment summaries
for af in sorted(Path("alignment_summaries").glob("*")):
    safe_copy(af, dirs["reports"] / "summaries" / af.name)

# Assembly selection method files
(dirs["reports"] / "summaries" / "assembly_selection").mkdir(exist_ok=True)
for mf in sorted(Path("method_files").glob("*")):
    safe_copy(mf, dirs["reports"] / "summaries" / "assembly_selection" / mf.name)

# QC - MultiQC
safe_copy("${multiqc_report}", dirs["reports"] / "qc" / "multiqc_report.html")

# Plasmid report
safe_copy("${plasmid_report}", dirs["reports"] / "summaries" / "plasmid_analysis.txt")

# =========================================================
# 03_results: Collated analysis results
# =========================================================
# Typing
safe_copy("${ompa_blast_tsv}", dirs["results"] / "typing" / "all_samples.ompa_blast.tsv")
safe_copy("${secondary_blast_tsv}", dirs["results"] / "typing" / "all_samples.secondary_blast.tsv")
safe_copy("${mlst_generic_tsv}", dirs["results"] / "typing" / "all_samples.mlst_chlamydiales.tsv")
safe_copy("${mlst_ct_tsv}", dirs["results"] / "typing" / "all_samples.mlst_ctrachomatis.tsv")

# Coverage
for cf in sorted(Path("coverage_files").glob("*")):
    safe_copy(cf, dirs["results"] / "coverage" / cf.name)

# Plasmid collated results
safe_copy("${plasmid_summary}", dirs["results"] / "plasmid" / "plasmid_summary.tsv")
safe_copy("${plasmid_blast_all}", dirs["results"] / "plasmid" / "plasmid_blast.tsv")
safe_copy("${plasmid_mlst_all}", dirs["results"] / "plasmid" / "plasmid_mlst.tsv")
safe_copy("${plasmid_report}", dirs["results"] / "plasmid" / "analysis_report.txt")

# =========================================================
# 04_phylogeny: Tree files
# =========================================================
tree_file = Path("${tree}")
if tree_file.exists() and tree_file.stat().st_size > 0:
    safe_copy(tree_file, dirs["phylogeny"] / f"{org}.tree")

# Tree data directory (SKA alignment and logs)
tree_data_dir = dirs["phylogeny"] / "tree_data"
tree_data_dir.mkdir(exist_ok=True)
# Copy any tree-related files that were staged
for tree_related in Path(".").glob("ska_alignment*"):
    safe_copy(tree_related, tree_data_dir / tree_related.name)
for tree_related in Path(".").glob("iqtree*"):
    safe_copy(tree_related, tree_data_dir / tree_related.name)

# =========================================================
# 05_annotation (if available)
# =========================================================
annotation_dir = Path("05_annotation")
annotation_dir.mkdir(exist_ok=True)
safe_copy("${annotation_summary}", annotation_dir / "annotation_summary.tsv")
safe_copy("${annotation_report}", annotation_dir / "annotation_report.txt")

# Coverage check TSV
safe_copy("${annotation_coverage_check}", annotation_dir / "coverage_check.tsv")

# Per-sample annotation symlinks
(annotation_dir / "per_sample").mkdir(exist_ok=True)
for sample_dir in sorted(dirs["per_sample"].iterdir()):
    if sample_dir.is_dir():
        sample_name = sample_dir.name
        # Check if annotation files exist in the main output
        # (they're published to {outdir}/{sample}/annotation/ by BAKTA_ANNOTATE)
        # We can't access those from here, but we note the structure

# =========================================================
# Mixed detection in reports/summaries
# =========================================================
safe_copy("${mixed_summary}", dirs["reports"] / "summaries" / "mixed_detection_summary.tsv")
safe_copy("${mixed_report}", dirs["reports"] / "summaries" / "mixed_detection_report.txt")

# =========================================================
# README
# =========================================================
# Count samples from assemblies
samples = sorted(set(
    d.name for d in dirs["per_sample"].iterdir() if d.is_dir()
))

readme_lines = [
    "# CtGAP-nf Output Directory",
    "Generated by CtGAP-nf Pipeline",
    "",
    "## Organized Structure",
    "",
    "01_per_sample/     - Sample-specific outputs",
    "   {sample}/",
    "       assembly.fasta     - Best assembly",
    "       *_ctgap_report.pdf - Per-sample report",
    "       *.plasmid.fasta    - Plasmid assembly",
    "",
    "02_reports/        - All reports in one place",
    "   per_sample/ - Individual sample PDFs",
    "   summaries/  - Cross-sample summary tables",
    "   qc/         - MultiQC report",
    "",
    "03_results/        - Collated analysis results",
    "   typing/     - ompA BLAST, MLST results (all samples)",
    "   coverage/   - Coverage statistics",
    "   plasmid/    - Plasmid BLAST, MLST, per-sample assemblies",
    "",
    "04_phylogeny/      - Phylogenetic tree (" + org + ".tree)",
    "",
    "05_annotation/     - Genome annotations (Bakta)",
    "   annotation_summary.tsv",
    "   annotation_report.txt",
    "",
    "## Key Files",
    "",
    "- Summary table:    02_reports/summaries/ctgap_all_samples.csv",
    "- Batch summary:    02_reports/summaries/batch_summary.txt",
    "- Genovar results:  03_results/typing/all_samples.ompa_blast.tsv",
    "- MLST results:     03_results/typing/all_samples.mlst_ctrachomatis.tsv",
    "- Mixed detection:  02_reports/summaries/mixed_detection_summary.tsv",
    "- Annotations:      05_annotation/annotation_summary.tsv",
    "- Phylogeny:        04_phylogeny/" + org + ".tree",
    "",
    "## Samples Processed: " + str(len(samples)),
    "",
]
for s in samples:
    readme_lines.append("  - " + s)
readme = "\\n".join(readme_lines) + "\\n"

with open("00_README.txt", "w") as f:
    f.write(readme)
    """
}
