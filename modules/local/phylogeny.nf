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

    # Log input list
    echo "=== SKA input list ===" > ska.log
    cat input.list >> ska.log
    echo "=== End input list ===" >> ska.log
    num_samples=\$(wc -l < input.list)
    echo "Number of input lines: \$num_samples" >> ska.log

    # Check if we have enough samples
    if [ "\$num_samples" -lt 2 ]; then
        echo "WARNING: Only \$num_samples samples in input list (need >=2 for SKA). Skipping." >> ska.log
        touch ska_alignment.fasta
    else
        set +e

        # Step 1: Build SKA sketch from input list
        echo "Running ska build..." >> ska.log
        ska build -o ska_out -k 17 -f input.list --threads ${task.cpus} >> ska.log 2>&1
        build_exit=\$?
        echo "ska build exit code: \$build_exit" >> ska.log

        if [ \$build_exit -ne 0 ] || [ ! -f ska_out.skf ]; then
            echo "WARNING: ska build failed. Creating empty alignment." >> ska.log
            touch ska_alignment.fasta
        else
            # Step 2: Map against reference with repeat masking
            echo "Running ska map..." >> ska.log
            ska map -o tmp_alignment --threads ${task.cpus} --repeat-mask ${params.scaffold_ref} ska_out.skf >> ska.log 2>&1
            map_exit=\$?
            echo "ska map exit code: \$map_exit" >> ska.log

            if [ \$map_exit -ne 0 ] || [ ! -f tmp_alignment ]; then
                echo "WARNING: ska map failed. Creating empty alignment." >> ska.log
                touch ska_alignment.fasta
            else
                # Step 3: Clean alignment (replace non-ACGTN with N, sanitize headers)
                echo "Cleaning alignment..." >> ska.log
                awk '
                    /^>/ {
                        # Clean header: remove special chars, keep first word
                        gsub(/#/, "_")
                        gsub(/:/, "_")
                        sub(/^>/, "")
                        split(\$0, a, " ")
                        print ">" a[1]
                        next
                    }
                    {
                        # Uppercase and replace non-ACGTN- with N
                        seq = toupper(\$0)
                        gsub(/[^ACGTN-]/, "N", seq)
                        print seq
                    }
                ' tmp_alignment > ska_alignment.fasta

                echo "Alignment created with \$(grep -c "^>" ska_alignment.fasta) sequences" >> ska.log

                # Cleanup
                rm -f ska_out.skf tmp_alignment
            fi
        fi

        set -e
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
