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
- **Pre-assembly coverage checks** — gates each requested assembly path at 10× mapped mean depth by default; routes low-coverage or empty-read samples to error reports
- **Host read removal** — HoCoRT + Bowtie2 depletion of human reads
- **Taxonomic filtering** — Kraken2 + KrakenTools extraction of Chlamydiales reads
- **Optional read subsampling** — seqtk caps extracted Chlamydiales read pairs per sample
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

Or use glob-based input. Both `{sample}_R{1,2}.fastq.gz` and Illumina-style
`{sample}_S##_L###_R{1,2}_001.fastq.gz` filenames are supported, including a mix
of the two styles in the same directory:

```bash
nextflow run main.nf --indir /path/to/fastqs/ --outdir results/
```

Read-mate suffixes (`_R1`/`_R2`, optionally followed by `_001`) and `.fastq.gz`
are always excluded from sample IDs. By default, `--truncate_sample_names true`
also removes trailing Illumina `_S##` or `_S##_L###` metadata:

| Example FASTQ | Default sample ID | With `--truncate_sample_names false` |
|---|---|---|
| `CHIR-018C_S84_L001_R1_001.fastq.gz` | `CHIR-018C` | `CHIR-018C_S84_L001` |
| `CHIR-030_S83_L001_R2_001.fastq.gz` | `CHIR-030` | `CHIR-030_S83_L001` |
| `V402426844_R1.fastq.gz` | `V402426844` | `V402426844` |

This applies only to glob-based input; CSV sample IDs are used as provided.
Each resulting sample ID must be unique. If multiple pairs collapse to the same
ID, the pipeline reports an error.

### 3. Assembly modes

```bash
# Auto mode (default) — recommended
nextflow run main.nf --input_csv samples.csv --mode auto

# De novo only
nextflow run main.nf --input_csv samples.csv --mode denovo

# Reference-guided + de novo
nextflow run main.nf --input_csv samples.csv --mode reference-denovo
```

### 4. Limit reads in high-depth samples

```bash
nextflow run main.nf --input_csv samples.csv -profile docker --subsample 1000000
```

`--subsample` sets the maximum number of **read pairs per sample** after
`EXTRACT_CHLAMYDIALES` (host removal and taxonomic filtering), before QC and assembly.
For example, `--subsample 1000000` retains up to 1,000,000 pairs (1,000,000 reads
in each of R1 and R2). Samples with fewer pairs retain all available pairs.
Provide a positive integer; fractions and zero are not accepted. Omit the option
to disable subsampling (the default).

The nf-core [`seqtk/sample`](https://nf-co.re/modules/seqtk_sample) module samples
both mates with the same fixed seed (`100`) to preserve pairing and make selection
reproducible. Two-pass sampling (`-2`) reduces memory usage and preserves read order.
QC and all downstream read-based analyses, including both assembly paths, plasmid
analysis, and mixed strain detection, use the sampled reads. These FASTQs are
published under `{outdir}/{sample}/subsample/`; the full extracted reads remain
under `{outdir}/{sample}/scrub/`.

---

## Parameters

All parameters can be set via command line (`--param value`) or in `nextflow.config`.

### Core Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_csv` | `null` | CSV samplesheet (columns: `sample`, `fastq_1`, `fastq_2`) |
| `--indir` | `../input` | Directory with `*_R{1,2}.fastq.gz` or `*_R{1,2}_001.fastq.gz` files (used if no CSV) |
| `--outdir` | `../output` | Output directory |
| `--mode` | `auto` | Assembly mode: `auto`, `denovo`, or `reference-denovo` |
| `--organism` | `ct` | Organism code |
| `--truncate_sample_names` | `true` | Remove trailing Illumina `_S##` or `_S##_L###` metadata from sample IDs; mate/file suffixes are always excluded (glob mode only) |

### Optional Features

| Parameter | Default | Description |
|---|---|---|
| `--annotation_enabled` | `true` | Run Bakta genome annotation |
| `--ct_typing_enabled` | `true` | Run Mash + FastANI strain typing |
| `--mixed_detection_enabled` | `true` | Run mixed strain detection |
| `--assembly_filter_enabled` | `true` | Filter non-Ct contigs from assemblies |
| `--subsample` | `null` (disabled) | Maximum extracted Chlamydiales read pairs per sample; positive integer, applied before QC and downstream analyses |

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
| `--min_assembly_depth` | `10` | Minimum mapped mean depth (×) required before each Shovill assembly path |
| `--ragtag_min_len` | `1000` | Minimum contig length for scaffolding |

See `nextflow.config` for the full list of parameters.

### Assembly Coverage Gate

```bash
nextflow run main.nf --input_csv samples.csv -profile docker --min_assembly_depth 10
```

Each requested assembly path must pass a coverage check before Shovill runs:

- **De novo:** after `DENOVO_SAMTOOLS_REF24`, use the highest per-reference mean
  depth from `samtools coverage` against the reference set. Unused references do
  not dilute the depth by being averaged together.
- **Reference-guided:** use coverage against the selected reference, calculated
  by `REF_SAMTOOLS_NAMESORT`. After `BAMTOFASTQ`, also check the actual extracted
  paired FASTQs destined for `REF_SHOVILL`, so an empty extraction is rejected
  even if the alignment contained mapped reads.

Depth must be at least `--min_assembly_depth` (default **10×**). Both mates must
contain reads and have equal read counts. Header-only coverage is treated as zero
depth. The threshold accepts finite, non-negative numbers; setting it to `0`
removes the positive-depth requirement but still rejects empty or unbalanced pairs.
Checks use post-subsampling reads when `--subsample` is enabled.
The Python coverage-check helper is staged into each task's work directory by
Nextflow, so coverage gates also run in remote AWS Batch containers.

In `auto` mode, both passing paths are scored normally. If only one passes, its
assembly is selected. In single-path modes, only that requested path is eligible.
Samples with no passing requested path skip assembly and assembly-dependent
analyses, including plasmid analysis, and receive error reports instead of a
normal genome report. Other samples continue. Low coverage is a sample-level
skip, while malformed FASTQs, invalid coverage files, and tool failures remain
task failures.

Coverage-check outputs are published directly, including when **all samples are
skipped** and no final assembly output is created:

- `01_per_sample/{sample}/{sample}.{path}.coverage_error.txt`: per-path error report.
- `02_reports/errors/{sample}.{path}.coverage_error.txt`: copies of error reports.
- `02_reports/summaries/assembly_coverage_checks.tsv`: all evaluated paths and metrics.
- `02_reports/summaries/assembly_errors.tsv`: rejected paths and reasons.
- `02_reports/summaries/assembly_coverage_summary.txt`: batch totals and samples
  with all requested paths skipped.
- `reports/coverage_checks/`: individual coverage-check TSVs.

Reports include mapped depth, the reference used, reference breadth, the threshold,
actual R1/R2 read counts, and reasons such as `NO_READ_PAIRS` or `LOW_MAPPED_DEPTH`.

---

## Output Structure

For samples that produce assemblies, the `REORGANIZE` step creates the output
structure below. Coverage-check summaries and per-path error reports are published
directly into the same numbered directories, even when all samples are skipped:

```
output/
├── 00_README.txt
├── 01_per_sample/
│   └── {sample}/
│       ├── assembly.fasta              # Best genome assembly
│       ├── {sample}_ctgap_report.pdf   # Comprehensive report
│       ├── {sample}_plasmid_report.pdf # Plasmid report
│       ├── {sample}.plasmid.fasta      # Plasmid assembly
│       └── {sample}.{path}.coverage_error.txt # If an assembly path was skipped
│
├── 02_reports/
│   ├── per_sample/                     # Individual PDFs
│   ├── errors/                         # Per-path coverage error reports
│   ├── summaries/
│   │   ├── ctgap_all_samples.csv       # Master summary
│   │   ├── batch_summary.txt
│   │   ├── assembly_methods.txt
│   │   ├── mixed_detection_summary.tsv
│   │   ├── assembly_coverage_checks.tsv # All evaluated paths and metrics
│   │   ├── assembly_errors.tsv          # Rejected paths and reasons
│   │   ├── assembly_coverage_summary.txt # Batch coverage-gate outcome
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

`{path}` is `denovo` or `reference-denovo`. A sample with no eligible requested
path receives coverage error reports instead of an assembly and normal genome
report. A sample with one passing path in `auto` mode can have both a successful
assembly and an error report for its rejected path.

---

## Assembly Selection Scoring (Auto Mode)

When using `--mode auto`, CtGAP scores each assembly on a 0–100 scale if both
paths pass the pre-assembly coverage gate. If only one path passes, its assembly
is selected without comparison; if neither passes, the sample is reported as skipped.

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
| 2e | `SEQTK_SAMPLE` | Optional subsampling of extracted read pairs (`--subsample`) |
| 3 | `FASTP_QC`, `MULTIQC` | Post-scrub QC reports |
| 4 | `DENOVO_BOWTIE_ALIGN` → `DENOVO_SAMTOOLS_REF24`, `SHOVILL`, `SCAFFOLD`, `GAPFILLER`, `RENAME_CONTIGS`, `QUAST` | De novo assembly path |
| 4b / 5e | `CHECK_DENOVO_COVERAGE`, `CHECK_REF_COVERAGE`, `ASSEMBLY_COVERAGE_SUMMARY` | Gate each requested assembly path on depth and read availability; report rejected paths |
| 5 | `REF_BOWTIE_ALIGN` → `REF_SAMTOOLS_NAMESORT`, `REF_BOWTIE_REF24_ALIGN` → `REF_SAMTOOLS_REF24`, `BAMTOFASTQ`, `REF_SHOVILL`, `REF_SCAFFOLD`, `REF_GAPFILLER`, `REF_RENAME`, `REF_QUAST` | Reference-guided assembly path |
| 6 | `SELECT_BEST_ASSEMBLY` / `COPY_BEST_ASSEMBLY` | Select the best eligible assembly; auto mode uses the sole passing path when the other is skipped |
| 6b | `FILTER_INDEX_ASSEMBLY` → `FILTER_MINIMAP2_ALIGN` → `FILTER_APPLY_CONTIGS` | Remove non-Ct contigs (optional) |
| 6c | `QUAST_FINAL` | QUAST on final (post-selection, post-filter) assembly for accurate report metrics |
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

### No reads remaining after host depletion

If Kraken2 completes successfully without producing classification records or a
report, the pipeline creates empty output files so `KRAKEN2_CLASSIFY` can emit its
required outputs. With no classification records, `EXTRACT_CHLAMYDIALES` skips
KrakenTools and produces valid empty gzip-compressed FASTQs. These steps therefore
handle samples with no remaining reads; the assembly coverage gates then route
those samples to error reports. A non-zero Kraken2 exit remains a task failure.

### Shovill fails or produces no contigs

`SHOVILL` and `REF_SHOVILL` preserve Shovill's exit status when logging through
`tee` and require a non-empty `shovill_tmp/contigs.fa` before publishing an assembly.
If a task fails, inspect `shovill.log` and `.command.out` in the task's Nextflow
work directory for the underlying Shovill/SPAdes error. Check that the input
FASTQs contain enough usable read pairs. The earlier error
`mv: can't rename 'shovill_tmp/contigs.fa': No such file or directory` can be a
secondary symptom of an assembly failure; the wrappers now stop at that failure
or report explicitly that no contigs were produced.

---

## Coverage Gate Checks

Run the Python metric checks and synthetic Nextflow routing tests from the
repository root (the routing tests do not run a real assembler or require databases):

```bash
python3 -B -m unittest discover -s tests -v
nextflow run tests/assembly_routing.nf --outdir coverage-test-output
nextflow run tests/assembly_routing.nf --all_rejected true --outdir coverage-test-empty-output
```
