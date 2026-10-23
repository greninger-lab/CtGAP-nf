/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_BLASTDB - Build BLAST databases for ompA and secondary typing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules: prepare_ompa_blastdb, prepare_secondary_blastdb
    (7-typing.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// PREPARE_OMPA_BLASTDB and PREPARE_SECONDARY_BLASTDB removed
// Pre-built BLAST databases used directly from mount point via
// params.ompa_blast_db_prefix and params.secondary_blast_db_prefix

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLAST_OMPA - ompA BLAST typing with optional secondary BLAST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: blast_ompa (7-typing.smk)
    Runs primary BLAST against ompA database. If top hit is A/B/Ba/C,
    runs secondary BLAST for finer resolution.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLAST_OMPA {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/blast:2.17.0--hb02a186_1"

    publishDir "${params.outdir}/${sample_id}/best/blast", mode: 'copy', pattern: '*.tab'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.blast.ompa.tab"),      emit: primary
    tuple val(sample_id), path("${sample_id}.blast.secondary.tab"),  emit: secondary

    script:
    """
    # Primary BLAST against ompA database (pre-built on mount)
    blastn \\
        -query ${assembly} \\
        -db ${params.ompa_blast_db_prefix} \\
        -num_threads ${task.cpus} \\
        -use_index false \\
        -max_target_seqs 1 \\
        -outfmt 6 \\
        -out ${sample_id}.blast.ompa.tab 2> blast.log || true

    # Check for empty results
    if [ ! -s ${sample_id}.blast.ompa.tab ]; then
        echo "WARNING: No BLAST hits found for ${sample_id}" >> blast.log
        touch ${sample_id}.blast.secondary.tab
        exit 0
    fi

    # Get top hit by bitscore (column 12)
    top_hit=\$(sort -k12,12nr ${sample_id}.blast.ompa.tab | head -n1 | cut -f2)

    # Check if top hit matches A/B/Ba/C genotypes requiring secondary BLAST
    target_hits="A_DQ064279.1 B_M33636.1 Ba_DQ064282.1 C_AF352789.1"
    if echo "\$target_hits" | grep -qw "\$top_hit"; then
        echo "Top hit \$top_hit matched - running secondary BLAST" >> blast.log
        blastn \\
            -query ${assembly} \\
            -db ${params.secondary_blast_db_prefix} \\
            -num_threads ${task.cpus} \\
            -use_index false \\
            -max_target_seqs 5 \\
            -outfmt 6 \\
            -out ${sample_id}.blast.secondary.tab 2>> blast.log || true
    else
        echo "Top hit \$top_hit - no secondary BLAST needed" >> blast.log
        touch ${sample_id}.blast.secondary.tab
    fi
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLATE_BLAST - Aggregate ompA BLAST results across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules: collate_blast, collate_secondary_blast (7-typing.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COLLATE_BLAST {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ('primary/*')
    path ('secondary/*')

    output:
    path "best.ompA_genovar.blast.tsv",  emit: primary_tsv
    path "best.secondary.blast.tsv",      emit: secondary_tsv

    script:
    """
    # Primary BLAST collation
    echo -e "query\\tsubject\\tpident\\tlength\\tmismatch\\tgapopen\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tevalue\\tbitscore" > best.ompA_genovar.blast.tsv
    for f in primary/*.tab; do [ -s "\$f" ] && cat "\$f" || true; done >> best.ompA_genovar.blast.tsv

    # Secondary BLAST collation
    echo -e "query\\tsubject\\tpident\\tlength\\tmismatch\\tgapopen\\tquery_start\\tquery_end\\tsubject_start\\tsubject_end\\tevalue\\tbitscore" > best.secondary.blast.tsv
    for f in secondary/*.tab; do [ -s "\$f" ] && cat "\$f" || true; done >> best.secondary.blast.tsv
    """
}
