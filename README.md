# CtGAP — Nextflow Port

***Chlamydia trachomatis* Genome Assembly Pipeline**

This is a **Nextflow DSL2 adaptation** of [CtGAP v2.0](https://github.com/D-Dean-Lab/CtGAP/tree/main), originally implemented as a Snakemake workflow by the [D-Dean-Lab](https://github.com/D-Dean-Lab). The pipeline has been adapted to Nextflow to run in our compute environment with Docker containers on AWS Batch.

> **Original repository:** <https://github.com/D-Dean-Lab/CtGAP>
>
> **Original authors:** Lucile Zhu, Ammar Aziz, Shola Olagoke, and contributors
>
> **License:** MIT (see [LICENSE](https://github.com/D-Dean-Lab/CtGAP/blob/main/LICENSE))

---

## Overview

CtGAP assembles *C. trachomatis* genomes from Illumina paired-end reads and performs:

- **Automatic assembly selection** — runs both de novo and reference-guided assembly, scores each, picks the best per sample
- **Host read removal** — HoCoRT + Bowtie2 depletion of human reads
- **Taxonomic filtering** — Kraken2 + KrakenTools extraction of Chlamydiales reads
- **Assembly** — Shovill/SPAdes + RagTag scaffolding + Gap2Seq gap filling
- **Assembly filtering** — removes non-Ct contigs via minimap2
- **ompA genotyping** — BLAST with secondary disambiguation (A/B/Ba/C)
- **MLST** — dual-scheme (Chlamydiales + *C. trachomatis*)
- **Ct strain typing** — whole-genome Mash + FastANI against reference genomes (optional)
- **Mixed strain detection** — heterozygous site analysis for co-infections (optional)
- **Dedicated plasmid module** — reference-guided extraction + Unicycler assembly + BLAST + MLST
- **Genome annotation** — Bakta annotation (optional)
- **Phylogenetics** — SKA2 alignment + IQ-TREE
- **Comprehensive PDF reports** — per-sample reports with all results
- **QC integration** — fastp, QUAST, MultiQC, coverage analysis

---

## Key Differences from the Snakemake Version

| Aspect | Snakemake (original) | Nextflow (this port) |
|---|---|---|
| **Workflow engine** | Snakemake ≥7.0 | Nextflow ≥23.04.0 (DSL2) |
| **Environment management** | Conda environments per rule | Docker containers per process |
| **Database setup** | `make setup` downloads databases | Databases pre-built on shared S3 mount point |
| **Index building** | Built at runtime, cached locally | Pre-built on mount point, no runtime index building |
| **Input format** | `{sample}_R{1,2}.fastq.gz` in `input/` dir | CSV samplesheet (`sample,fastq_1,fastq_2`) or glob-based |
| **Execution** | `make run` / `./ctgap run` | `nextflow run main.nf` |
| **Compute environment** | Local / HPC with conda | AWS Batch with Docker containers |
| **Process isolation** | Single tool per conda env | Single tool per Docker container (bowtie2, samtools, etc. run in separate containers) |

---

## Requirements

### Software

- **Nextflow** ≥23.04.0
- **Docker** (each process runs in its own container — no tools need to be installed on `PATH`)

All bioinformatics tools are provided via Docker containers specified in each process definition:

| Category | Tools | Container |
|---|---|---|
| **Read QC** | fastp | `quay.io/biocontainers/fastp` |
| **Read QC** | MultiQC | `quay.io/biocontainers/multiqc` |
| **Host removal** | HoCoRT | `quay.io/biocontainers/hocort` |
| **Host removal** | Bowtie2 | `quay.io/biocontainers/bowtie2` |
| **Host removal** | SAMtools | `quay.io/biocontainers/samtools` |
| **Classification** | Kraken2 | `quay.io/biocontainers/kraken2` |
| **Classification** | KrakenTools | `quay.io/jefffurlong/krakentools` |
| **Assembly** | Shovill (SPAdes) | `quay.io/biocontainers/shovill` |
| **Assembly** | Unicycler | `quay.io/biocontainers/unicycler` |
| **Assembly** | RagTag | `quay.io/biocontainers/ragtag` |
| **Assembly** | Gap2Seq | `quay.io/biocontainers/gap2seq` |
| **Alignment** | minimap2 | `quay.io/biocontainers/minimap2` |
| **Alignment** | bedtools | `quay.io/biocontainers/bedtools` |
| **Typing** | BLAST+ | `quay.io/biocontainers/blast` |
| **Typing** | claMLST (pyMLST) | `quay.io/biocontainers/pymlst` |
| **Strain typing** | Mash | `quay.io/biocontainers/mash` |
| **Strain typing** | FastANI | `quay.io/biocontainers/fastani` |
| **Annotation** | Bakta | `quay.io/biocontainers/bakta` |
| **Variant calling** | BCFtools | `quay.io/biocontainers/bcftools` |
| **Sequence tools** | seqkit | `quay.io/biocontainers/seqkit` |
| **Sequence tools** | seqtk | `quay.io/biocontainers/seqtk` |
| **Phylogenetics** | Augur (IQ-TREE) | `quay.io/biocontainers/augur` |
| **Reporting** | Python 3 + pandas + reportlab | `quay.io/jefffurlong/pandas-reportlab` |

### Databases (pre-provisioned on shared mount point)

All databases and reference indices are **pre-built** and available on a shared S3 mount point accessible from AWS Batch containers. No database building occurs at runtime.

| Database | Default Path | Purpose |
|---|---|---|
| GRCh38 reference | `resources/grch38.fasta` | Host read depletion |
| GRCh38 bt2 index | `resources/grch38.bt2/grch38.*.bt2` | Pre-built bowtie2 index for host depletion |
| Kraken2 standard | `resources/standardDB/` | Taxonomic classification |
| Bakta light DB | `resources/bakta_db/db-light/` | Genome annotation |
| Ct references | `resources/references/ct/` | ompA, MLST, refset, scaffold, individual genomes, plasmid refs |
| Refset bt2 index | `resources/references/ct/refset.*.bt2` | Pre-built bowtie2 index for coverage alignment |
| Plurality bt2 index | `resources/references/ct/plurality.*.bt2` | Pre-built bowtie2 index for ref-guided alignment |
| Plasmid bt2 index | `resources/references/ct/pCT.*.bt2` | Pre-built bowtie2 index for plasmid read mapping |
| ompA BLAST DB | `resources/references/ct/ompaBlastDB/ompa21.*` | Pre-built BLAST database for ompA typing |
| Secondary BLAST DB | `resources/references/ct/secondaryBlastDB/secondary.*` | Pre-built BLAST database for A/B/Ba/C disambiguation |
| Plasmid BLAST DB | `resources/references/ct/20_Ct_plasmids.*` | Pre-built BLAST database for plasmid typing |
| Mash sketch | `resources/references/ct/ct_refs.msh` | Pre-built Mash sketch for CT strain typing |
| MLST databases | `resources/references/ct/pubmlst/` | claMLST schemes (Chlamydiales, C. trachomatis, plasmid) |

Input FASTQ files (specified via `--input_csv` or `--indir`) are staged by Nextflow automatically — they can be local paths, relative paths, or S3 URIs. Only the databases and reference genomes listed above must be accessible on the shared mount point inside AWS Batch containers.


---

## Quick Start

### 1. Prepare a samplesheet

Create a CSV file with a header row:

```csv
sample,fastq_1,fastq_2
CT001,/path/to/CT001_R1.fastq.gz,/path/to/CT001_R2.fastq.gz
CT002,/path/to/CT002_R1.fastq.gz,/path/to/CT002_R2.fastq.gz
```

### 2. Run the pipeline

```bash
nextflow run main.nf --input_csv samples.csv --outdir results/
```

FASTQ paths in the CSV can be local paths, relative paths, or S3 URIs — Nextflow stages them automatically for each job (including AWS Batch). Only the reference databases and genomes need to be on the shared mount point.

Or use glob-based input (files must be named `{sample}_R{1,2}.fastq.gz`):

```bash
nextflow run main.nf --indir /path/to/fastqs/ --outdir results/
```

> **Note:** By default, Illumina `_S##_L###` suffixes are stripped from sample names
> (e.g. `CT_Cap_1_S31_L001` becomes `CT_Cap_1`). Use `--truncate_sample_names false`
> to keep the original names.

### 3. Assembly modes

```bash
# Auto mode (default) — recommended
nextflow run main.nf --input_csv samples.csv --mode auto

# De novo only
nextflow run main.nf --input_csv samples.csv --mode denovo

# Reference-guided + de novo
nextflow run main.nf --input_csv samples.csv --mode reference-denovo
```

---

## Parameters

All parameters can be set via command line (`--param value`) or in `nextflow.config`.

### Core Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_csv` | `null` | CSV samplesheet (columns: `sample`, `fastq_1`, `fastq_2`) |
| `--indir` | `../input` | Directory with `*_R{1,2}.fastq.gz` files (used if no CSV) |
| `--outdir` | `../output` | Output directory |
| `--mode` | `auto` | Assembly mode: `auto`, `denovo`, or `reference-denovo` |
| `--organism` | `ct` | Organism code |
| `--truncate_sample_names` | `true` | Strip Illumina `_S##_L###` suffixes from sample names (glob mode only) |

### Optional Features

| Parameter | Default | Description |
|---|---|---|
| `--annotation_enabled` | `true` | Run Bakta genome annotation |
| `--ct_typing_enabled` | `false` | Run Mash + FastANI strain typing |
| `--mixed_detection_enabled` | `true` | Run mixed strain detection |
| `--assembly_filter_enabled` | `true` | Filter non-Ct contigs from assemblies |

### Database Paths

| Parameter | Default | Description |
|---|---|---|
| `--host_reference` | `/mnt/.../grch38.fasta` | GRCh38 human reference |
| `--kraken_db` | `/mnt/.../standardDB` | Kraken2 database |
| `--bakta_db` | `/mnt/.../bakta_db` | Bakta annotation database |
| `--resources_path` | `/mnt/.../references` | Ct reference genomes directory |

### Assembly Parameters

| Parameter | Default | Description |
|---|---|---|
| `--shovill_gsize` | `1.04m` | Expected genome size |
| `--shovill_depth` | `200` | Shovill depth limit |
| `--ragtag_min_len` | `1000` | Minimum contig length for scaffolding |

See `nextflow.config` for the full list of parameters.

---

## Output Structure

After the pipeline completes, the `REORGANIZE` step creates a clean output structure:

```
output/
├── 00_README.txt
├── 01_per_sample/
│   └── {sample}/
│       ├── assembly.fasta              # Best genome assembly
│       ├── {sample}_ctgap_report.pdf   # Comprehensive report
│       ├── {sample}_plasmid_report.pdf # Plasmid report
│       └── {sample}.plasmid.fasta      # Plasmid assembly
│
├── 02_reports/
│   ├── per_sample/                     # Individual PDFs
│   ├── summaries/
│   │   ├── ctgap_all_samples.csv       # Master summary
│   │   ├── batch_summary.txt
│   │   ├── assembly_methods.txt
│   │   ├── mixed_detection_summary.tsv
│   │   └── assembly_selection/         # Per-sample method files
│   └── qc/
│       └── multiqc_report.html
│
├── 03_results/
│   ├── typing/
│   │   ├── all_samples.ompa_blast.tsv
│   │   └── all_samples.mlst_ctrachomatis.tsv
│   ├── coverage/
│   └── plasmid/
│       ├── plasmid_summary.tsv
│       └── per_sample/
│
├── 04_phylogeny/
│   └── ct.tree                         # Newick format
│
└── 05_annotation/                      # If annotation enabled
    ├── annotation_summary.tsv
    └── annotation_report.txt
```

---

## Assembly Selection Scoring (Auto Mode)

When using `--mode auto`, CtGAP scores each assembly on a 0–100 scale:

| Metric | Weight | Full Score Criteria |
|---|---|---|
| N50 | 40 pts | ≥95% of expected genome size (1.04 Mb) |
| Contig count | 30 pts | ≤2 contigs |
| Length accuracy | 20 pts | Within 10 kb of expected 1.04 Mb |
| GC content | 10 pts | 40–42% (optimal for Ct) |

De novo assembly is preferred in ties or close calls (within 2 points), as it introduces less reference bias.

---

## Pipeline Steps

| Step | Process(es) | Description |
|---|---|---|
| 1 | `FASTP` | Adapter trimming and quality filtering |
| 2 | `HOCORT_DEPLETE`, `KRAKEN2_CLASSIFY`, `EXTRACT_CHLAMYDIALES` | Host removal + taxonomic filtering |
| 3 | `FASTP_QC`, `MULTIQC` | Post-scrub QC reports |
| 4 | `DENOVO_BOWTIE_ALIGN` → `DENOVO_SAMTOOLS_REF24`, `SHOVILL`, `SCAFFOLD`, `GAPFILLER`, `RENAME_CONTIGS`, `QUAST` | De novo assembly path |
| 5 | `REF_BOWTIE_ALIGN` → `REF_SAMTOOLS_NAMESORT`, `REF_BOWTIE_REF24_ALIGN` → `REF_SAMTOOLS_REF24`, `BAMTOFASTQ`, `REF_SHOVILL`, `REF_SCAFFOLD`, `REF_GAPFILLER`, `REF_RENAME`, `REF_QUAST` | Reference-guided assembly path |
| 6 | `SELECT_BEST_ASSEMBLY` / `COPY_BEST_ASSEMBLY` | Assembly selection (auto mode scores both, picks best) |
| 6b | `FILTER_INDEX_ASSEMBLY` → `FILTER_MINIMAP2_ALIGN` → `FILTER_APPLY_CONTIGS` | Remove non-Ct contigs (optional) |
| 7 | `BLAST_OMPA`, `MLST`, `SKA_ALIGNMENT`, `IQTREE` | Typing + phylogenetics |
| 8 | `CT_MASH_DIST`, `CT_FASTANI`, `CT_ASSIGN` | Ct strain typing (optional) |
| 9 | `PLASMID_BOWTIE_ALIGN` → `PLASMID_SAMTOOLS_SORT`, `EXTRACT_PLASMID_READS`, `UNICYCLER_PLASMID`, `ANALYZE_PLASMID`, `REORIENT_PLASMID`, `BLAST_PLASMID`, `MLST_PLASMID`, `PLASMID_REPORT` | Dedicated plasmid assembly |
| 10 | `CHECK_ANNOTATION_COVERAGE`, `BAKTA_ANNOTATE`, `ANNOTATION_SUMMARY` | Genome annotation (optional, coverage-gated) |
| 11 | `MIXED_EXTRACT_MAIN_CONTIG`, `MIXED_BOWTIE_INDEX` → `MIXED_BOWTIE_ALIGN` → `MIXED_SAMTOOLS_SORT` → `MIXED_BCFTOOLS_CALL` | Mixed strain detection (optional) |
| 12 | `COMPREHENSIVE_REPORT`, `SUMMARY_CSV`, `BATCH_SUMMARY` | Reporting |
| 13 | `REORGANIZE` | Clean output structure |

> **Note:** Arrows (→) indicate processes that are split for Docker container isolation. For example, `DENOVO_BOWTIE_ALIGN` (bowtie2 container) produces a SAM file that is passed to `DENOVO_SAMTOOLS_REF24` (samtools container) for sorting, indexing, and coverage calculation.

---

## Profiles

```bash
# Run with Docker (local)
nextflow run main.nf --input_csv samples.csv -profile docker

# Run on AWS Batch with Docker
nextflow run main.nf --input_csv samples.csv -profile docker_batch

# Run with Singularity
nextflow run main.nf --input_csv samples.csv -profile singularity

# Local testing (overrides database paths to local copies)
nextflow run main.nf --input_csv samples.csv -profile local

# Test profile (uses test data paths)
nextflow run main.nf -profile test
```

The `docker_batch` profile is the primary execution mode, running on AWS Batch with Docker containers and databases on an S3 mount point. The `local` profile overrides database paths for local development/testing.

---

## Troubleshooting

### "No valid samples found"

Check your CSV format — it must have a header row with columns `sample`, `fastq_1`, `fastq_2`:

```csv
sample,fastq_1,fastq_2
mysample,/full/path/to/reads_R1.fastq.gz,/full/path/to/reads_R2.fastq.gz
```

### Pipeline failed mid-run

Nextflow automatically caches completed tasks. Just re-run with `-resume`:

```bash
nextflow run main.nf --input_csv samples.csv -resume
```

### Database not found

Ensure mount points are accessible inside the AWS Batch containers. The mount point base path is configured in `nextflow.config` as `params.resources_base`. Verify the path:

```bash
ls -la /mnt/s3/dlmp-research/labs/greninger/greninger-lab-file-share/refs/CtGAP/resources/
```

If running locally, use the `-profile local` to override database paths to local copies.

### Low assembly quality

Check the per-sample PDF report for coverage statistics. Low input depth (<10x) affects assembly contiguity.

---

