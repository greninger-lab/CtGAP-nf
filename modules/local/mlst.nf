/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MLST - claMLST typing (Chlamydiales + C. trachomatis schemes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: mlst (7-typing.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MLST {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pymlst:2.3.2--pyhdfd78af_0"

    publishDir "${params.outdir}/${sample_id}/best/mlst", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.genome.chlamydiales.mlst.txt"),  emit: generic
    tuple val(sample_id), path("${sample_id}.genome.ctrachomatis.mlst.txt"),  emit: ct

    script:
    """
    claMLST search \\
        ${params.mlst_db}/chlamydiales \\
        ${assembly} > ${sample_id}.genome.chlamydiales.mlst.txt 2> mlst.log

    claMLST search \\
        ${params.mlst_db}/c.trachomatis \\
        ${assembly} > ${sample_id}.genome.ctrachomatis.mlst.txt 2>> mlst.log
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLATE_MLST - Aggregate MLST results across all samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: collate_mlst (7-typing.smk)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process COLLATE_MLST {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ('generic/*')
    path ('ct/*')

    output:
    path "best.mlst.generic.results.tsv", emit: generic
    path "best.mlst.ct.results.tsv",      emit: ct

    script:
    """
    # Collate generic MLST
    first_generic=\$(ls generic/* | head -1)
    head -1 "\$first_generic" > best.mlst.generic.results.tsv
    for f in generic/*; do
        tail -n +2 "\$f" >> best.mlst.generic.results.tsv
    done

    # Collate CT MLST
    first_ct=\$(ls ct/* | head -1)
    head -1 "\$first_ct" > best.mlst.ct.results.tsv
    for f in ct/*; do
        tail -n +2 "\$f" >> best.mlst.ct.results.tsv
    done
    """
}
