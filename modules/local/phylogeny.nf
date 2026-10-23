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

    container "quay.io/biocontainers/ska2:0.5.1--h4349ce8_0"

    publishDir "${params.outdir}/tree", mode: 'copy'

    input:
    path ('assemblies/*')

    output:
    path "ska_alignment.fasta", emit: alignment
    path "ska.log",             emit: log

    script:
    """
    # Build input list: sample_name<tab>fasta_path
    > input.list
    for fasta in assemblies/*; do
        sample=\$(basename "\$fasta" | sed 's/[.]assembly[.]filtered[.]fasta\$//' | sed 's/[.]assembly[.]fasta\$//' | sed 's/[.]fasta\$//' | sed 's/[.]fa\$//')
        echo -e "\$sample\\t\$fasta" >> input.list
    done

    # Add individual reference fastas from mount point if present
    if [ -d "${params.individual_refs}" ] && [ "\$(ls -A ${params.individual_refs} 2>/dev/null)" ]; then
        for fasta in ${params.individual_refs}/*; do
            sample=\$(basename "\$fasta" | sed 's/[.]fasta\$//' | sed 's/[.]fa\$//')
            echo -e "\$sample\\t\$fasta" >> input.list
        done
    fi

    # Debug: show input list and environment
    echo "=== SKA input list ===" >> ska.log
    cat input.list >> ska.log
    echo "=== End input list ===" >> ska.log
    echo "Number of input lines: \$(wc -l < input.list)" >> ska.log
    echo "Number of assembly files: \$(ls assemblies/ | wc -l)" >> ska.log
    echo "Assembly files:" >> ska.log
    ls -la assemblies/ >> ska.log 2>&1
    echo "which ska: \$(which ska 2>&1)" >> ska.log
    echo "which python3: \$(which python3 2>&1)" >> ska.log
    echo "which generate_ska_alignment.py: \$(which generate_ska_alignment.py 2>&1)" >> ska.log
    echo "which gubbins_alignment_checker.py: \$(which gubbins_alignment_checker.py 2>&1)" >> ska.log
    echo "PATH: \$PATH" >> ska.log
    echo "=== End debug ===" >> ska.log

    # Check if we have enough samples
    num_samples=\$(wc -l < input.list)
    if [ "\$num_samples" -lt 2 ]; then
        echo "WARNING: Only \$num_samples samples in input list (need >=2 for SKA). Skipping." >> ska.log
        touch ska_alignment.fasta
    else
        # Run SKA alignment
        set +e
        generate_ska_alignment.py \\
            --input input.list \\
            --out ska_alignment.fasta \\
            --reference ${params.scaffold_ref} \\
            --threads ${task.cpus} 2>&1 | tee -a ska.log
        exitcode=\$?
        set -e

        # If SKA fails, create empty alignment
        if [ \$exitcode -ne 0 ] || [ ! -f ska_alignment.fasta ]; then
            echo "WARNING: SKA alignment failed (exit code \$exitcode). Creating empty alignment." >> ska.log
            cat ska.log >&2
            touch ska_alignment.fasta
        fi
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
