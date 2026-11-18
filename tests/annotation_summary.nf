nextflow.enable.dsl = 2
params.outdir = 'annotation-summary-test-output'
params.annotated = false
params.annotation_min_coverage = 90
params.ct_reference_size = 1000

include { ANNOTATION_SUMMARY; ANNOTATION_COVERAGE_CHECK } from '../modules/local/annotate'

process GENERATE {
    output:
    path 'sample.quast.report.tsv', emit: quast
    path 'sample.txt', emit: annotation, optional: true
    script:
    def n50 = params.annotated ? 1000 : 500
    """
    printf 'Assembly\\ttest\\nN50\\t${n50}\\n' > sample.quast.report.tsv
    if [ '${params.annotated}' = 'true' ]; then
        printf 'Sequences: 1\\nSize: 1000\\nGC: 41.0\\nCDSs: 5\\ntRNAs: 1\\nrRNAs: 0\\nncRNAs: 0\\nhypotheticals: 1\\npseudogenes: 0\\n' > sample.txt
    fi
    """
}

workflow {
    GENERATE()
    ANNOTATION_COVERAGE_CHECK(GENERATE.out.quast.collect())
    def no_file = file("${projectDir}/../assets/NO_FILE", checkIfExists: true)
    annotations = GENERATE.out.annotation.collect().ifEmpty { no_file }
    ANNOTATION_SUMMARY(annotations, ANNOTATION_COVERAGE_CHECK.out.coverage_check)
    ANNOTATION_SUMMARY.out.summary.map { tsv ->
        def rows = tsv.text.readLines()
        assert rows.size() == (params.annotated ? 2 : 1)
        if (params.annotated) {
            assert rows[1] == 'sample\t1\t1000\t41.0\t5\t1\t0\t0\t1\t0'
        }
        "PASS: annotation TSV contains ${params.annotated ? 1 : 0} real samples"
    }.view()
    ANNOTATION_SUMMARY.out.report.map { report ->
        assert report.text.contains("Samples annotated: ${params.annotated ? 1 : 0}")
        assert report.text.contains(params.annotated ? 'sample\t1000\t100.0\tYES' : 'sample\t500\t50.0\tNO')
        if (!params.annotated) {
            assert report.text.contains('No samples met the coverage threshold for annotation.')
        }
        'PASS: annotation report preserves coverage checks and sample counts'
    }.view()
}
