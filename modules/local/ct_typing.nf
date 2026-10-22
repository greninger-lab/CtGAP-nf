/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CT_TYPING - ANI + Mash-based CT reference typing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Equivalent to Snakemake rules from 6-ct_typing.smk:
      ct_sketch_refs, ct_mash_sketch_query, ct_mash_dist, ct_mash_best,
      ct_fastani, ct_fastani_best, ct_assign, ct_typing_summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// CT_SKETCH_REFS removed — pre-built Mash sketch used directly
// from mount point via params.ct_refs_mash_sketch

/*
 * Per-sample: Mash sketch + distance + best hit
 */
process CT_MASH_DIST {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/mash:2.3--hf85e966_11"
    
    publishDir "${params.outdir}/${sample_id}/ct_typing/mash", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_vs_refs.dist.tsv"),  emit: dist
    tuple val(sample_id), path("${sample_id}_best_mash.tsv"),     emit: best

    script:
    def k = params.ct_mash_k
    def s = params.ct_mash_sketch_size
    """
    # Sketch query
    mash sketch -o ${sample_id} -k ${k} -s ${s} ${assembly} 2> sketch.log

    # Calculate distances (using pre-built sketch from mount point)
    mash dist ${params.ct_refs_mash_sketch} ${sample_id}.msh | sort -k3,3g > ${sample_id}_vs_refs.dist.tsv 2> dist.log

    # Select best hit (already sorted by mash_dist ascending)
    echo -e "ref\\tquery\\tmash_dist\\tpvalue\\tshared_hashes\\tsample" > ${sample_id}_best_mash.tsv
    head -1 ${sample_id}_vs_refs.dist.tsv | awk -F'\\t' -v sample="${sample_id}" 'BEGIN{OFS="\\t"} {print \$1, \$2, \$3, \$4, \$5, sample}' >> ${sample_id}_best_mash.tsv
    """
}

/*
 * Per-sample: fastANI against all CT references
 */
process CT_FASTANI {
    tag "${sample_id}"
    label 'process_medium'

    container "quay.io/biocontainers/fastani:1.34--hb66fcc3_7"
    
    publishDir "${params.outdir}/${sample_id}/ct_typing/ani", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_vs_refs.ani.tsv"), emit: ani

    script:
    def fragLen = params.ct_fastani_fragLen
    def minFraction = params.ct_fastani_minFraction
    """
    # Create file lists from mount point
    ls ${params.individual_refs}/*.fa* > ref_list.txt
    echo "${assembly}" > query_list.txt

    # Run fastANI
    fastANI --rl ref_list.txt --ql query_list.txt \\
            -o ${sample_id}_vs_refs.ani.tsv \\
            --fragLen ${fragLen} \\
            --minFraction ${minFraction} \\
            --threads ${task.cpus} 2> fastani.log || touch ${sample_id}_vs_refs.ani.tsv
    """
}

/*
 * Per-sample: Final CT type assignment (combines ANI + Mash)
 */
process CT_ASSIGN {
    tag "${sample_id}"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/${sample_id}/ct_typing", mode: 'copy'

    input:
    tuple val(sample_id), path(ani_tsv), path(mash_best), path(mash_dist)

    output:
    tuple val(sample_id), path("${sample_id}.ct_assignment.csv"), emit: assignment

    script:
    """
    #!/usr/bin/env python3
    import os
    import pandas as pd
    from pathlib import Path

    sample = "${sample_id}"
    delta_ani = ${params.ct_delta_ani}
    min_cov = ${params.ct_min_ani_cov}
    mash_margin = ${params.ct_min_mash_margin}
    high_ani = ${params.ct_high_ani}
    med_ani = ${params.ct_med_ani}

    # Read fastANI results
    ani_file = "${ani_tsv}"
    if os.path.getsize(ani_file) == 0:
        result = pd.DataFrame([{
            "sample": sample, "ref": "NO_ANI_RESULT", "ANI": pd.NA,
            "frags_over": 0, "frags_total": 0, "cov": 0.0,
            "alt_ref": "", "alt_ANI": pd.NA, "alt_cov": pd.NA,
            "flag": "no_ani_result"
        }])
    else:
        df = pd.read_csv(ani_file, sep="\\t", header=None,
                         names=["query", "ref", "ANI", "frags_over", "frags_total"])
        df["ANI"] = pd.to_numeric(df["ANI"], errors='coerce')
        df["frags_over"] = pd.to_numeric(df["frags_over"], errors='coerce')
        df["frags_total"] = pd.to_numeric(df["frags_total"], errors='coerce')
        df["cov"] = df["frags_over"] / df["frags_total"]
        df = df.sort_values(["ANI", "cov"], ascending=[False, False]).reset_index(drop=True)

        top = df.iloc[0]
        result = pd.DataFrame([{
            "sample": sample, "ref": top["ref"], "ANI": top["ANI"],
            "frags_over": top["frags_over"], "frags_total": top["frags_total"],
            "cov": top["cov"], "alt_ref": "", "alt_ANI": pd.NA,
            "alt_cov": pd.NA, "flag": "clear"
        }])

        if top["cov"] < min_cov:
            result.loc[0, "flag"] = "low_ani_coverage"

        if len(df) > 1:
            alt = df.iloc[1]
            result.loc[0, "alt_ref"] = alt["ref"]
            result.loc[0, "alt_ANI"] = alt["ANI"]
            result.loc[0, "alt_cov"] = alt["cov"]

            if (top["ANI"] - alt["ANI"]) < delta_ani:
                result.loc[0, "flag"] = "near_tie"
                # Use Mash to break tie
                mash_df = pd.read_csv("${mash_dist}", sep="\\t", header=None,
                                      names=["ref", "query", "mash_dist", "pvalue", "shared_hashes"])
                mash_df["mash_dist"] = pd.to_numeric(mash_df["mash_dist"], errors='coerce')
                mash_df["ref_base"] = mash_df["ref"].apply(lambda x: os.path.basename(str(x)))
                top_ref = os.path.basename(str(top["ref"]))
                alt_ref = os.path.basename(str(alt["ref"]))
                m_top = mash_df[mash_df["ref_base"] == top_ref]["mash_dist"].min()
                m_alt = mash_df[mash_df["ref_base"] == alt_ref]["mash_dist"].min()
                if pd.notna(m_top) and pd.notna(m_alt):
                    if (m_alt - m_top) > mash_margin:
                        result.loc[0, "flag"] = "clear_by_mash_vote"
                    elif (m_top - m_alt) > mash_margin:
                        result.loc[0, "ref"] = alt["ref"]
                        result.loc[0, "ANI"] = alt["ANI"]
                        result.loc[0, "cov"] = alt["cov"]
                        result.loc[0, "flag"] = "alt_by_mash_vote"

    # Read Mash best for merging
    mash = pd.read_csv("${mash_best}", sep="\\t")
    result = result.merge(mash.rename(columns={"ref": "ref_mash"}), on="sample", how="left")

    def clean_ref(ref_path):
        if pd.isna(ref_path) or ref_path == "NO_ANI_RESULT":
            return ref_path
        name = Path(str(ref_path)).name
        for ext in [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"]:
            if name.endswith(ext):
                return name[:-len(ext)]
        return Path(str(ref_path)).stem

    result["ref_clean"] = result["ref"].apply(clean_ref)
    result["ref_mash_clean"] = result["ref_mash"].apply(clean_ref)

    def get_confidence(row):
        ani_val = row["ANI"]
        cov_val = row["cov"]
        flag = row["flag"]
        if pd.isna(ani_val) or flag == "no_ani_result":
            return "Low/Needs review (no ANI)"
        if cov_val < min_cov:
            return "Low/Needs review (low coverage)"
        if flag == "near_tie":
            return "Review (near-tie)"
        if ani_val >= high_ani:
            return f"High (ANI>={high_ani}%)"
        elif ani_val >= med_ani:
            return f"Medium ({med_ani}-{high_ani}%)"
        return "Low/Needs review"

    result["confidence"] = result.apply(get_confidence, axis=1)
    result["putative_reference"] = result.apply(
        lambda r: r["ref_clean"] if (pd.notna(r["ANI"]) and r["cov"] >= min_cov)
                  else r["ref_mash_clean"], axis=1)

    result.to_csv("${sample_id}.ct_assignment.csv", index=False)
    """
}

/*
 * Aggregate CT typing results for all samples
 */
process CT_TYPING_SUMMARY {
    tag "all_samples"
    label 'process_low'

    container "quay.io/biocontainers/pandas:2.2.1"

    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    path ('assignments/*')

    output:
    path "ct_typing_all_samples.csv", emit: summary

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    from pathlib import Path

    dfs = []
    for f in sorted(Path("assignments").glob("*.csv")):
        dfs.append(pd.read_csv(f))

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.sort_values("sample")
    combined.to_csv("ct_typing_all_samples.csv", index=False)
    """
}
