#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CtGAP - Chlamydia trachomatis Genome Assembly Pipeline (Nextflow)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ported from the Snakemake implementation.
    Uses DSL2 with modular process definitions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print pipeline header
log.info """\
    ===============================================
       C t G A P - n f
       Chlamydia trachomatis Genome Assembly (Nextflow)
    ===============================================
    mode         : ${params.mode}
    input        : ${params.input_csv ?: params.indir}
    output       : ${params.outdir}
    organism     : ${params.organism}
    reference    : ${params.reference}
    host ref     : ${params.host_reference}
    truncate IDs : ${params.truncate_sample_names}
    ===============================================
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INCLUDE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                } from './modules/local/fastp'
include { HOCORT_DEPLETE       } from './modules/local/hocort_deplete'
include { KRAKEN2_CLASSIFY     } from './modules/local/kraken2'
include { EXTRACT_CHLAMYDIALES } from './modules/local/extract_chlamydiales'
include { FASTP_QC             } from './modules/local/fastp_qc'
include { MULTIQC              } from './modules/local/multiqc'
include { DENOVO_BOWTIE_ALIGN    } from './modules/local/denovo_bowtie_ref24'
include { DENOVO_SAMTOOLS_REF24  } from './modules/local/denovo_bowtie_ref24'
include { SHOVILL              } from './modules/local/shovill'
include { SCAFFOLD             } from './modules/local/scaffold'
include { GAPFILLER            } from './modules/local/gapfiller'
include { RENAME_CONTIGS       } from './modules/local/rename_contigs'
include { QUAST                } from './modules/local/quast'
include { DENOVO_COLLATE_COVERAGE } from './modules/local/denovo_collate'
include { DENOVO_COLLATE_FLAGSTAT } from './modules/local/denovo_collate'
include { SCAFFOLD_SUMMARY        } from './modules/local/denovo_collate'
include { REF_BOWTIE_ALIGN         } from './modules/local/ref_bowtie'
include { REF_SAMTOOLS_SORT as REF_SAMTOOLS_NAMESORT } from './modules/local/ref_bowtie'
include { REF_BOWTIE_REF24_ALIGN   } from './modules/local/ref_bowtie_ref24'
include { REF_SAMTOOLS_REF24       } from './modules/local/ref_bowtie_ref24'
include { BAMTOFASTQ               } from './modules/local/bamtofastq'
include { REF_SHOVILL              } from './modules/local/ref_shovill'
include { REF_SCAFFOLD             } from './modules/local/ref_scaffold'
include { REF_GAPFILLER            } from './modules/local/ref_gapfiller'
include { REF_RENAME               } from './modules/local/ref_rename'
include { REF_QUAST                } from './modules/local/ref_quast'
include { REF_QUAST_SCAFFOLDED     } from './modules/local/ref_quast'
include { REF_COLLATE_COVERAGE     } from './modules/local/ref_collate'
include { REF_COLLATE_FLAGSTAT     } from './modules/local/ref_collate'
include { REF_SCAFFOLD_SUMMARY     } from './modules/local/ref_collate'
include { SELECT_BEST_ASSEMBLY     } from './modules/local/select_best'
include { COPY_BEST_ASSEMBLY       } from './modules/local/select_best'
include { ASSEMBLY_METHOD_SUMMARY  } from './modules/local/select_best'
include { FILTER_INDEX_ASSEMBLY    } from './modules/local/filter_contigs'
include { FILTER_MINIMAP2_ALIGN    } from './modules/local/filter_contigs'
include { FILTER_APPLY_CONTIGS     } from './modules/local/filter_contigs'
include { FILTER_SUMMARY           } from './modules/local/filter_contigs'
include { BLAST_OMPA               } from './modules/local/blast_typing'
include { COLLATE_BLAST            } from './modules/local/blast_typing'
include { MLST                     } from './modules/local/mlst'
include { COLLATE_MLST             } from './modules/local/mlst'
include { SKA_ALIGNMENT            } from './modules/local/phylogeny'
include { IQTREE                   } from './modules/local/phylogeny'
include { CT_MASH_DIST             } from './modules/local/ct_typing'
include { CT_FASTANI               } from './modules/local/ct_typing'
include { CT_ASSIGN                } from './modules/local/ct_typing'
include { CT_TYPING_SUMMARY        } from './modules/local/ct_typing'
include { PLASMID_BOWTIE_ALIGN     } from './modules/local/plasmid'
include { PLASMID_SAMTOOLS_SORT    } from './modules/local/plasmid'
include { EXTRACT_PLASMID_READS    } from './modules/local/plasmid'
include { UNICYCLER_PLASMID        } from './modules/local/plasmid'
include { ANALYZE_PLASMID          } from './modules/local/plasmid'
include { REORIENT_PLASMID         } from './modules/local/plasmid'
include { BLAST_PLASMID            } from './modules/local/plasmid'
include { MLST_PLASMID             } from './modules/local/plasmid'
include { COLLATE_PLASMID          } from './modules/local/plasmid'
include { PLASMID_REPORT           } from './modules/local/plasmid'
include { CHECK_ANNOTATION_COVERAGE } from './modules/local/annotate'
include { BAKTA_ANNOTATE           } from './modules/local/annotate'
include { ANNOTATION_SUMMARY       } from './modules/local/annotate'
include { ANNOTATION_COVERAGE_CHECK } from './modules/local/annotate'
include { MIXED_EXTRACT_MAIN_CONTIG } from './modules/local/mixed_detection'
include { MIXED_BOWTIE_INDEX       } from './modules/local/mixed_detection'
include { MIXED_BOWTIE_ALIGN       } from './modules/local/mixed_detection'
include { MIXED_SAMTOOLS_SORT      } from './modules/local/mixed_detection'
include { MIXED_BCFTOOLS_CALL      } from './modules/local/mixed_detection'
include { MIXED_DETECTION_SUMMARY  } from './modules/local/mixed_detection'
include { COMPREHENSIVE_REPORT     } from './modules/local/reports'
include { SUMMARY_CSV              } from './modules/local/reports'
include { BATCH_SUMMARY            } from './modules/local/reports'
include { REORGANIZE               } from './modules/local/reorganize'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // ---- Build input channel from paired FASTQ files ----
    if (params.input_csv) {
        // CSV input with header row: sample,fastq_1,fastq_2
        // fastq_2 is optional -- if missing, inferred as single-end
        // This pipeline requires paired-end reads
        ch_reads = Channel
            .fromPath(params.input_csv, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                def sample_id = row.sample?.trim()
                if (!sample_id) {
                    error "Missing 'sample' column or empty sample name in CSV: ${row}"
                }
                def r1_path = row.fastq_1?.trim()
                def r2_path = row.fastq_2?.trim()
                if (!r1_path) {
                    error "Missing 'fastq_1' for sample '${sample_id}' in ${params.input_csv}"
                }
                if (!r2_path) {
                    error "Single-end reads not supported. Please provide 'fastq_2' for sample '${sample_id}'"
                }
                def r1 = file(r1_path)
                def r2 = file(r2_path)
                tuple(sample_id, r1, r2)
            }
            .ifEmpty { error "No valid samples found in ${params.input_csv}" }
            // Check for duplicate sample IDs
            .groupTuple()
            .map { sample_id, r1_list, r2_list ->
                if (r1_list.size() > 1) {
                    error "Duplicate sample '${sample_id}' found in ${params.input_csv}. Each sample must have a unique name."
                }
                tuple(sample_id, r1_list[0], r2_list[0])
            }
    } else {
        // Glob-based input: expects {sample}_R{1,2}.fastq.gz in params.indir
        ch_reads = Channel
            .fromFilePairs("${params.indir}/*_R{1,2}.fastq.gz", size: 2)
            .ifEmpty { error "No paired FASTQ files found in ${params.indir} matching *_R{1,2}.fastq.gz" }
            .map { sample_id, reads ->
                if (reads.size() != 2) {
                    error "Expected paired-end files (R1 and R2) for sample '${sample_id}', but found ${reads.size()} file(s)"
                }
                // Optionally strip Illumina _S##_L### suffixes from sample names
                def clean_id = sample_id
                if (params.truncate_sample_names) {
                    clean_id = sample_id.replaceAll(/_S\d+(_L\d+)?$/, '')
                    if (clean_id != sample_id) {
                        log.info "  Sample name truncated: ${sample_id} -> ${clean_id}"
                    }
                }
                tuple(clean_id, reads[0], reads[1])
            }
            // Check for duplicate sample names (can occur after truncation)
            .groupTuple()
            .map { sample_id, r1_list, r2_list ->
                if (r1_list.size() > 1) {
                    error "Duplicate sample name '${sample_id}' found in ${params.indir} " +
                          "(${r1_list.size()} file pairs). " +
                          (params.truncate_sample_names
                              ? "This may be caused by --truncate_sample_names stripping _S##_L### suffixes. " +
                                "Use --truncate_sample_names false to keep original names, or rename input files."
                              : "Each sample must have a unique name.")
                }
                tuple(sample_id, r1_list[0], r2_list[0])
            }
    }

    // ---- Reference files ----
    // Note: Mount-point databases/references are accessed via params.* directly
    // inside process scripts, not as Nextflow path channels, because mount points
    // are only available inside AWS Batch containers at runtime.

    // ---- STEP 1: Trim reads with fastp ----
    FASTP(ch_reads)

    // ---- STEP 2: Host read scrubbing ----
    // 2a. Deplete host reads (uses pre-built bt2 index from mount point)
    HOCORT_DEPLETE(FASTP.out.reads)

    // 2c. Classify remaining reads with Kraken2
    KRAKEN2_CLASSIFY(HOCORT_DEPLETE.out.reads)

    // 2d. Extract Chlamydiales reads
    // Combine host-depleted reads with kraken2 output for extraction
    ch_for_extract = HOCORT_DEPLETE.out.reads
        .join(KRAKEN2_CLASSIFY.out.results)
        .map { sample_id, r1, r2, k2_output, k2_report ->
            tuple(sample_id, r1, r2, k2_output, k2_report)
        }
    EXTRACT_CHLAMYDIALES(ch_for_extract)

    // Expose scrubbed reads for downstream steps
    ch_scrubbed = EXTRACT_CHLAMYDIALES.out.reads

    // ---- STEP 3: QC on scrubbed reads ----
    // 3a. Per-sample fastp QC (report only, no trimming)
    FASTP_QC(ch_scrubbed)

    // 3b. Aggregate all QC reports with MultiQC
    MULTIQC(FASTP_QC.out.json.collect())

    // ---- STEP 4: De novo assembly ----
    // 4a. Align scrubbed reads to reference set for coverage stats (pre-built bt2 index)
    DENOVO_BOWTIE_ALIGN(ch_scrubbed)
    DENOVO_SAMTOOLS_REF24(DENOVO_BOWTIE_ALIGN.out.sam)

    // 4c. De novo assembly with Shovill
    SHOVILL(ch_scrubbed)

    // 4d. Scaffold contigs with RagTag
    SCAFFOLD(SHOVILL.out.contigs)

    // 4e. Fill gaps with Gap2Seq
    // Combine scaffold output with scrubbed reads for gap filling
    ch_for_gapfill = SCAFFOLD.out.scaffold
        .join(ch_scrubbed)
        .map { sample_id, scaffold, r1, r2 ->
            tuple(sample_id, scaffold, r1, r2)
        }
    GAPFILLER(ch_for_gapfill)

    // 4f. Rename contigs with sample-prefixed names
    RENAME_CONTIGS(GAPFILLER.out.filled)

    // 4g. Assembly quality statistics
    QUAST(RENAME_CONTIGS.out.assembly)

    // 4h. Collate coverage and flagstat across all samples
    ch_all_coverages = DENOVO_SAMTOOLS_REF24.out.coverage
        .map { sample_id, cov -> cov }
        .collect()
    DENOVO_COLLATE_COVERAGE(ch_all_coverages)

    ch_all_flagstats = DENOVO_SAMTOOLS_REF24.out.flagstat
        .map { sample_id, fs -> fs }
        .collect()
    DENOVO_COLLATE_FLAGSTAT(ch_all_flagstats)

    // 4i. Scaffold summary
    SCAFFOLD_SUMMARY(SCAFFOLD.out.result.collect())

    // Expose de novo assembly for downstream steps
    ch_denovo_assembly = RENAME_CONTIGS.out.assembly

    // ---- STEP 5: Reference-guided assembly (runs in 'reference-denovo' and 'auto' modes) ----
    if (params.mode == 'reference-denovo' || params.mode == 'auto') {

        // 5a. Align scrubbed reads to specified reference (pre-built bt2 index, name-sorted for bamtofastq)
        REF_BOWTIE_ALIGN(ch_scrubbed)
        REF_SAMTOOLS_NAMESORT(REF_BOWTIE_ALIGN.out.sam)

        // 5b. Align scrubbed reads to reference set for coverage stats (pre-built bt2 index)
        REF_BOWTIE_REF24_ALIGN(ch_scrubbed)
        REF_SAMTOOLS_REF24(REF_BOWTIE_REF24_ALIGN.out.sam)

        // 5d. Extract aligned reads as FASTQ
        BAMTOFASTQ(REF_SAMTOOLS_NAMESORT.out.bam)

        // 5e. Assemble extracted reads with Shovill
        REF_SHOVILL(BAMTOFASTQ.out.reads)

        // 5f. Scaffold with RagTag
        REF_SCAFFOLD(REF_SHOVILL.out.contigs)

        // 5g. Gap filling with Gap2Seq (uses scrubbed reads)
        ch_ref_for_gapfill = REF_SCAFFOLD.out.scaffold
            .join(ch_scrubbed)
            .map { sample_id, scaffold, r1, r2 ->
                tuple(sample_id, scaffold, r1, r2)
            }
        REF_GAPFILLER(ch_ref_for_gapfill)

        // 5h. Rename contigs
        REF_RENAME(REF_GAPFILLER.out.filled)

        // 5i. Assembly quality statistics (reference-free)
        REF_QUAST(REF_RENAME.out.assembly)

        // 5j. Assembly quality statistics vs scaffold reference
        REF_QUAST_SCAFFOLDED(REF_RENAME.out.assembly)

        // 5k. Collate coverage and flagstat
        ch_ref_coverages = REF_SAMTOOLS_REF24.out.coverage
            .map { sample_id, cov -> cov }
            .collect()
        REF_COLLATE_COVERAGE(ch_ref_coverages)

        ch_ref_flagstats = REF_SAMTOOLS_NAMESORT.out.flagstat
            .map { sample_id, fs -> fs }
            .collect()
        REF_COLLATE_FLAGSTAT(ch_ref_flagstats)

        // 5l. Scaffold summary
        REF_SCAFFOLD_SUMMARY(REF_SCAFFOLD.out.result.collect())

        ch_refguided_assembly = REF_RENAME.out.assembly
    }

    // ---- STEP 6: Select best assembly ----
    if (params.mode == 'auto') {
        // Join denovo and ref-guided assemblies + their QUAST reports by sample_id
        ch_for_select = ch_denovo_assembly
            .join(ch_refguided_assembly)
            .join(QUAST.out.report)
            .join(REF_QUAST.out.report)
            .map { sample_id, denovo_asm, refguided_asm, denovo_quast, refguided_quast ->
                tuple(sample_id, denovo_asm, refguided_asm, denovo_quast, refguided_quast)
            }
        SELECT_BEST_ASSEMBLY(ch_for_select)
        ch_best_assembly = SELECT_BEST_ASSEMBLY.out.assembly
        ch_method_files  = SELECT_BEST_ASSEMBLY.out.method

    } else if (params.mode == 'denovo') {
        COPY_BEST_ASSEMBLY(ch_denovo_assembly, 'denovo')
        ch_best_assembly = COPY_BEST_ASSEMBLY.out.assembly
        ch_method_files  = COPY_BEST_ASSEMBLY.out.method

    } else if (params.mode == 'reference-denovo') {
        COPY_BEST_ASSEMBLY(ch_refguided_assembly, 'reference-denovo')
        ch_best_assembly = COPY_BEST_ASSEMBLY.out.assembly
        ch_method_files  = COPY_BEST_ASSEMBLY.out.method

    } else {
        error "Invalid mode '${params.mode}'. Must be: auto, denovo, or reference-denovo"
    }

    // 6b. Assembly method summary
    ASSEMBLY_METHOD_SUMMARY(ch_method_files.collect(), params.mode)

    // ---- STEP 6b: Assembly filtering (optional) ----
    if (params.assembly_filter_enabled) {
        FILTER_INDEX_ASSEMBLY(ch_best_assembly)
        FILTER_MINIMAP2_ALIGN(FILTER_INDEX_ASSEMBLY.out.indexed)
        FILTER_APPLY_CONTIGS(FILTER_MINIMAP2_ALIGN.out.aligned)

        // Use filtered assembly for downstream
        ch_final_assembly = FILTER_APPLY_CONTIGS.out.assembly

        // Collate filter reports
        FILTER_SUMMARY(FILTER_APPLY_CONTIGS.out.report.collect())
    } else {
        ch_final_assembly = ch_best_assembly
    }

    // ---- STEP 7: Typing ----
    // 7a. BLAST ompA typing (uses pre-built BLAST DBs from mount point)
    BLAST_OMPA(ch_final_assembly)

    // 7c. MLST typing
    MLST(ch_final_assembly)

    // 7d. Collate BLAST results
    ch_primary_tabs = BLAST_OMPA.out.primary
        .map { sample_id, tab -> tab }
        .collect()
    ch_secondary_tabs = BLAST_OMPA.out.secondary
        .map { sample_id, tab -> tab }
        .collect()
    COLLATE_BLAST(ch_primary_tabs, ch_secondary_tabs)

    // 7e. Collate MLST results
    ch_mlst_generic = MLST.out.generic
        .map { sample_id, txt -> txt }
        .collect()
    ch_mlst_ct = MLST.out.ct
        .map { sample_id, txt -> txt }
        .collect()
    COLLATE_MLST(ch_mlst_generic, ch_mlst_ct)

    // 7f. Phylogenetic tree
    ch_assemblies_for_tree = ch_final_assembly
        .map { sample_id, fasta -> fasta }
        .collect()

    SKA_ALIGNMENT(ch_assemblies_for_tree)
    IQTREE(SKA_ALIGNMENT.out.alignment)

    // ---- STEP 8: CT Typing (optional) ----
    if (params.ct_typing_enabled) {
        // Per-sample Mash distance (uses pre-built sketch from mount)
        CT_MASH_DIST(ch_final_assembly)

        // Per-sample fastANI
        CT_FASTANI(ch_final_assembly)

        // Per-sample assignment (combine ANI + Mash)
        ch_for_assign = CT_FASTANI.out.ani
            .join(CT_MASH_DIST.out.best)
            .join(CT_MASH_DIST.out.dist)
        CT_ASSIGN(ch_for_assign)

        // Aggregate all assignments
        ch_all_assignments = CT_ASSIGN.out.assignment
            .map { sample_id, csv -> csv }
            .collect()
        CT_TYPING_SUMMARY(ch_all_assignments)
    }

    // ---- STEP 9: Dedicated Plasmid Assembly ----
    // 9a. Align scrubbed reads to reference plasmid (uses pre-built bt2 index from mount)
    PLASMID_BOWTIE_ALIGN(ch_scrubbed)
    PLASMID_SAMTOOLS_SORT(PLASMID_BOWTIE_ALIGN.out.sam)

    // 9c. Extract plasmid read pairs
    EXTRACT_PLASMID_READS(PLASMID_SAMTOOLS_SORT.out.bam)

    // 9d. Unicycler assembly
    ch_for_unicycler = EXTRACT_PLASMID_READS.out.reads
        .join(EXTRACT_PLASMID_READS.out.count)
        .map { sample_id, r1, r2, count_file ->
            tuple(sample_id, r1, r2, count_file)
        }
    UNICYCLER_PLASMID(ch_for_unicycler)

    // 9e. Analyze plasmid assembly
    ANALYZE_PLASMID(
        UNICYCLER_PLASMID.out.assembly,
        PLASMID_SAMTOOLS_SORT.out.stats,
        EXTRACT_PLASMID_READS.out.count
    )

    // 9f. Re-orient plasmid
    REORIENT_PLASMID(ANALYZE_PLASMID.out.plasmid)

    // 9g. BLAST plasmid (uses pre-built BLAST DB from mount point)
    BLAST_PLASMID(REORIENT_PLASMID.out.reoriented)

    // 9h. MLST plasmid
    MLST_PLASMID(ANALYZE_PLASMID.out.plasmid)

    // 9i. Collate plasmid results
    ch_plasmid_analyses = ANALYZE_PLASMID.out.analysis
        .map { sample_id, txt -> txt }
        .collect()
    ch_plasmid_blasts = BLAST_PLASMID.out.blast
        .map { sample_id, tsv -> tsv }
        .collect()
    ch_plasmid_mlsts = MLST_PLASMID.out.mlst
        .map { sample_id, txt -> txt }
        .collect()
    COLLATE_PLASMID(ch_plasmid_analyses, ch_plasmid_blasts, ch_plasmid_mlsts)

    // 9j. Per-sample plasmid PDF reports
    ch_for_plasmid_report = ANALYZE_PLASMID.out.analysis
        .join(BLAST_PLASMID.out.blast)
        .join(MLST_PLASMID.out.mlst)
        .join(PLASMID_SAMTOOLS_SORT.out.stats)
    PLASMID_REPORT(ch_for_plasmid_report)

    // ---- STEP 10: Genome Annotation (optional, coverage-gated) ----
    if (params.annotation_enabled) {
        // Check coverage for each sample using QUAST reports
        CHECK_ANNOTATION_COVERAGE(QUAST.out.report)

        // Filter to only samples that pass coverage threshold, then join with assemblies
        ch_passing = CHECK_ANNOTATION_COVERAGE.out.result
            .filter { sample_id, passes -> passes == "true" }
            .map { sample_id, passes -> tuple(sample_id) }

        ch_for_bakta = ch_final_assembly
            .join(ch_passing)

        BAKTA_ANNOTATE(ch_for_bakta)

        // Generate coverage check TSV for all samples
        ch_quast_for_cov = QUAST.out.report
            .map { sample_id, report -> report }
            .collect()
        ANNOTATION_COVERAGE_CHECK(ch_quast_for_cov)

        // Collate annotation results
        ch_annotation_txts = BAKTA_ANNOTATE.out.txt
            .map { sample_id, txt -> txt }
            .collect()
            .ifEmpty { file('NO_ANNOTATIONS') }

        ANNOTATION_SUMMARY(ch_annotation_txts, ANNOTATION_COVERAGE_CHECK.out.coverage_check)
    }

    // ---- STEP 11: Mixed Strain Detection (optional) ----
    if (params.mixed_detection_enabled) {
        // Extract main contig from each assembly
        MIXED_EXTRACT_MAIN_CONTIG(ch_final_assembly)

        // Combine main contig with scrubbed reads for variant calling
        ch_for_mixed = MIXED_EXTRACT_MAIN_CONTIG.out.main_contig
            .join(ch_scrubbed)
            .map { sample_id, main_contig, r1, r2 ->
                tuple(sample_id, main_contig, r1, r2)
            }

        // Split into separate container steps
        MIXED_BOWTIE_INDEX(ch_for_mixed)
        MIXED_BOWTIE_ALIGN(MIXED_BOWTIE_INDEX.out.indexed)
        MIXED_SAMTOOLS_SORT(MIXED_BOWTIE_ALIGN.out.aligned)
        MIXED_BCFTOOLS_CALL(MIXED_SAMTOOLS_SORT.out.sorted)

        // Collate results
        ch_mixed_summaries = MIXED_BCFTOOLS_CALL.out.summary
            .map { sample_id, txt -> txt }
            .collect()
        ch_mixed_het = MIXED_BCFTOOLS_CALL.out.het_sites
            .map { sample_id, tsv -> tsv }
            .collect()

        MIXED_DETECTION_SUMMARY(ch_mixed_summaries, ch_mixed_het)
    }

    // ---- STEP 12: Comprehensive Reporting ----
    // Build per-sample annotation channel (optional)
    if (params.annotation_enabled) {
        // Bakta .txt output tupled by sample_id; samples that didn't pass get a dummy
        ch_annotation_for_report = BAKTA_ANNOTATE.out.txt
            .map { sample_id, txt -> tuple(sample_id, txt) }
        // Right-join with all samples so missing ones get a placeholder
        ch_annotation_per_sample = ch_final_assembly
            .map { sample_id, fasta -> tuple(sample_id) }
            .join(ch_annotation_for_report, remainder: true)
            .map { items ->
                def sample_id = items[0]
                def txt = items.size() > 1 && items[1] != null ? items[1] : file('NO_ANNOTATION')
                tuple(sample_id, txt)
            }
    } else {
        ch_annotation_per_sample = ch_final_assembly
            .map { sample_id, fasta -> tuple(sample_id, file('NO_ANNOTATION')) }
    }

    // Build per-sample mixed detection channel (optional)
    if (params.mixed_detection_enabled) {
        ch_mixed_for_report = MIXED_BCFTOOLS_CALL.out.summary
    } else {
        ch_mixed_for_report = ch_final_assembly
            .map { sample_id, fasta -> tuple(sample_id, file('NO_MIXED')) }
    }

    // Build per-sample CT typing channel (optional)
    if (params.ct_typing_enabled) {
        ch_ct_typing_for_report = CT_ASSIGN.out.assignment
    } else {
        ch_ct_typing_for_report = ch_final_assembly
            .map { sample_id, fasta -> tuple(sample_id, file('NO_CT_TYPING')) }
    }

    // Build per-sample report input by joining key outputs
    ch_report_base = ch_final_assembly
        .join(QUAST.out.report)
        .join(BLAST_OMPA.out.primary)
        .join(BLAST_OMPA.out.secondary)
        .join(MLST.out.generic)
        .join(MLST.out.ct)
        .join(ANALYZE_PLASMID.out.analysis)
        .join(BLAST_PLASMID.out.blast)
        .join(MLST_PLASMID.out.mlst)
        .join(FASTP_QC.out.reports.map { sample_id, html, json -> tuple(sample_id, json) })
        .join(KRAKEN2_CLASSIFY.out.results.map { sample_id, k2, report -> tuple(sample_id, report) })
        .join(DENOVO_SAMTOOLS_REF24.out.coverage)
        .join(ch_annotation_per_sample)
        .join(ch_mixed_for_report)
        .join(ch_ct_typing_for_report)

    // Add method file by extracting sample_id from filename
    // Method files are named {sample_id}.method.txt
    ch_method_tupled = ch_method_files.map { f ->
        def name = f.name.toString().replace('.method.txt', '')
        tuple(name, f)
    }

    ch_report_input = ch_method_tupled
        .join(ch_report_base)
        .map { sample_id, method_file, assembly, quast, blast_p, blast_s,
               mlst_gen, mlst_ct, plasmid_analysis, plasmid_blast, plasmid_mlst,
               fastp_json, kraken_report, coverage, annotation_txt, mixed_summary,
               ct_typing_csv ->
            tuple(sample_id, assembly, method_file, quast, blast_p, blast_s,
                  mlst_gen, mlst_ct, plasmid_analysis, plasmid_blast, plasmid_mlst,
                  fastp_json, kraken_report, coverage, annotation_txt, mixed_summary,
                  ct_typing_csv)
        }

    COMPREHENSIVE_REPORT(ch_report_input)

    // ---- STEP 13: Summary CSV and Batch Summary ----
    // Collect inputs for summary CSV
    ch_all_methods = ch_method_files.collect()
    ch_all_quast = QUAST.out.report
        .map { sample_id, report -> report }
        .collect()
    ch_all_blast_tabs = BLAST_OMPA.out.primary
        .map { sample_id, tab -> tab }
        .collect()
    ch_all_mlst_ct = MLST.out.ct
        .map { sample_id, txt -> txt }
        .collect()
    ch_all_plasmid_analyses = ANALYZE_PLASMID.out.analysis
        .map { sample_id, txt -> txt }
        .collect()

    SUMMARY_CSV(
        ch_all_methods,
        ch_all_quast,
        ch_all_blast_tabs,
        ch_all_mlst_ct,
        ch_all_plasmid_analyses
    )

    BATCH_SUMMARY(SUMMARY_CSV.out.csv)

    // ---- STEP 14: Reorganize output ----
    ch_all_pdfs = COMPREHENSIVE_REPORT.out.pdf
        .map { sample_id, pdf -> pdf }
        .collect()
    ch_all_assemblies = ch_final_assembly
        .map { sample_id, fasta -> fasta }
        .collect()
    ch_all_plasmid_fastas = ANALYZE_PLASMID.out.plasmid
        .map { sample_id, fasta -> fasta }
        .collect()
    ch_all_plasmid_pdfs = PLASMID_REPORT.out.pdf
        .map { sample_id, pdf -> pdf }
        .collect()

    // Collect coverage files (denovo + optionally ref-guided)
    ch_coverage_for_reorg = DENOVO_COLLATE_COVERAGE.out.coverage
        .mix(
            params.mode != 'denovo'
                ? REF_COLLATE_COVERAGE.out.coverage
                : Channel.empty()
        )
        .collect()
        .ifEmpty { file('NO_COVERAGE') }

    // Collect scaffold summaries
    ch_scaffold_summaries = SCAFFOLD_SUMMARY.out.report
        .mix(
            params.mode != 'denovo'
                ? REF_SCAFFOLD_SUMMARY.out.report
                : Channel.empty()
        )
        .collect()
        .ifEmpty { file('NO_SCAFFOLD_SUMMARIES') }

    // Collect alignment summaries
    ch_alignment_summaries = DENOVO_COLLATE_FLAGSTAT.out.summary
        .mix(
            params.mode != 'denovo'
                ? REF_COLLATE_FLAGSTAT.out.summary
                : Channel.empty()
        )
        .collect()
        .ifEmpty { file('NO_ALIGNMENT_SUMMARIES') }

    // Annotation outputs (conditional)
    if (params.annotation_enabled) {
        ch_annotation_summary_file = ANNOTATION_SUMMARY.out.summary
        ch_annotation_report_file  = ANNOTATION_SUMMARY.out.report
        ch_annotation_coverage_file = ANNOTATION_COVERAGE_CHECK.out.coverage_check
    } else {
        ch_annotation_summary_file = Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: false).ifEmpty(file('NO_ANNOTATION_SUMMARY'))
        ch_annotation_report_file  = Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: false).ifEmpty(file('NO_ANNOTATION_REPORT'))
        ch_annotation_coverage_file = Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: false).ifEmpty(file('NO_ANNOTATION_COVERAGE'))
    }

    // Mixed detection outputs (conditional)
    if (params.mixed_detection_enabled) {
        ch_mixed_summary_file = MIXED_DETECTION_SUMMARY.out.summary_tsv
        ch_mixed_report_file  = MIXED_DETECTION_SUMMARY.out.report
    } else {
        ch_mixed_summary_file = Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: false).ifEmpty(file('NO_MIXED_SUMMARY'))
        ch_mixed_report_file  = Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: false).ifEmpty(file('NO_MIXED_REPORT'))
    }

    REORGANIZE(
        SUMMARY_CSV.out.csv,
        BATCH_SUMMARY.out.report,
        ch_all_pdfs,
        ch_all_assemblies,
        IQTREE.out.tree,
        ch_all_plasmid_fastas,
        ch_all_plasmid_pdfs,
        COLLATE_PLASMID.out.summary,
        COLLATE_PLASMID.out.blast_all,
        COLLATE_PLASMID.out.mlst_all,
        COLLATE_PLASMID.out.report,
        COLLATE_BLAST.out.primary_tsv,
        COLLATE_BLAST.out.secondary_tsv,
        COLLATE_MLST.out.generic,
        COLLATE_MLST.out.ct,
        ch_coverage_for_reorg,
        ASSEMBLY_METHOD_SUMMARY.out.summary,
        ch_method_files.collect(),
        MULTIQC.out.report,
        ch_scaffold_summaries,
        ch_alignment_summaries,
        ch_annotation_summary_file,
        ch_annotation_report_file,
        ch_annotation_coverage_file,
        ch_mixed_summary_file,
        ch_mixed_report_file
    )
}
