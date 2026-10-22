/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SKA_ALIGNMENT - Generate SKA alignment from assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules: ska_input_prep + ska_alignment (7-typing.smk)
    Prepares input list and runs SKA alignment against reference.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SKA_ALIGNMENT {
    tag "all_samples"
    label 'process_medium'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/tree", mode: 'copy'

    input:
    path ('assemblies/*')

    output:
    path "ska_alignment.fasta", emit: alignment

    script:
    """
    # Build input list: sample_name<tab>fasta_path
    > input.list
    for fasta in assemblies/*; do
        sample=\$(basename "\$fasta" | sed 's/\\.fasta\$//' | sed 's/\\.fa\$//')
        echo -e "\$sample\\t\$fasta" >> input.list
    done

    # Add individual reference fastas from mount point if present
    if [ -d "${params.individual_refs}" ] && [ "\$(ls -A ${params.individual_refs} 2>/dev/null)" ]; then
        for fasta in ${params.individual_refs}/*; do
            sample=\$(basename "\$fasta" | sed 's/\\.fasta\$//' | sed 's/\\.fa\$//')
            echo -e "\$sample\\t\$fasta" >> input.list
        done
    fi

    # Run SKA alignment
    set +e
    generate_ska_alignment.py \\
        --input input.list \\
        --out ska_alignment.fasta \\
        --reference ${params.scaffold_ref} \\
        --threads ${task.cpus} 2>&1 | tee ska.log
    exitcode=\$?
    set -e

    # If SKA fails (e.g., too few samples), create empty alignment
    if [ \$exitcode -ne 0 ] || [ ! -f ska_alignment.fasta ]; then
        echo "WARNING: SKA alignment failed (exit code \$exitcode). Creating empty alignment." >> ska.log
        touch ska_alignment.fasta
    fi
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IQTREE - Build phylogenetic tree from SKA alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rule: tree (7-typing.smk)
    Uses augur tree (IQ-TREE backend) with bootstrap support.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process IQTREE {
    tag "phylogeny"
    label 'process_medium'

    container "quay.io/biocontainers/augur:34.2.1--pyhdfd78af_0"

    publishDir "${params.outdir}", mode: 'copy', pattern: '*.tree'
    publishDir "${params.outdir}/tree", mode: 'copy', pattern: '*.log'

    input:
    path alignment

    output:
    path "${params.organism}.tree", emit: tree
    path "iqtree.log",              emit: log

    script:
    """
    # Check if alignment is non-empty and has enough sequences
    seq_count=\$(grep "^>" ${alignment} 2>/dev/null | wc -l | tr -d ' ')
    seq_count=\${seq_count:-0}

    if [ "\$seq_count" -lt 3 ]; then
        echo "WARNING: Alignment has \$seq_count sequences (need >=3 for tree). Skipping." > iqtree.log
        echo "(no_tree)" > ${params.organism}.tree
    else
        augur tree \\
            --alignment ${alignment} \\
            --output ${params.organism}.tree \\
            --override-default-args \\
            --tree-builder-args="-alrt 1000 -B 1000" \\
            --nthreads ${task.cpus} > iqtree.log 2>&1
    fi
    """
}
