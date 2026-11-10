nextflow.enable.dsl = 2
params.min_assembly_depth = 10
params.outdir = 'coverage-test-output'
params.all_rejected = false

include { CHECK_ASSEMBLY_COVERAGE as CHECK_DENOVO_COVERAGE; CHECK_ASSEMBLY_COVERAGE as CHECK_REF_COVERAGE; ASSEMBLY_COVERAGE_SUMMARY } from '../modules/local/assembly_coverage'
include { SELECT_BEST_ASSEMBLY } from '../modules/local/select_best'
include { CREATE_CANDIDATE as DENOVO_CANDIDATE; CREATE_CANDIDATE as REF_CANDIDATE } from './fixtures/assembly_candidate'

process GENERATE {
    input:
    val sample_id
    output:
    tuple val(sample_id), path('d_R1.fastq.gz'), path('d_R2.fastq.gz'), path('d.coverage.tsv'), emit: denovo
    tuple val(sample_id), path('r_R1.fastq.gz'), path('r_R2.fastq.gz'), path('r.coverage.tsv'), emit: reference
    script:
    """
    python3 - <<'PY'
    import gzip
    sample = '${sample_id}'
    for method in ('d', 'r'):
        for mate in (1, 2):
            with gzip.open(f'{method}_R{mate}.fastq.gz', 'wt') as handle:
                if not (sample == 'empty_ref' and method == 'r'):
                    handle.write('@read\\nACGT\\n+\\nIIII\\n')
        depth = 1 if sample == 'none' or (sample == 'denovo_only' and method == 'r') or (sample == 'ref_only' and method == 'd') else 20
        with open(f'{method}.coverage.tsv', 'w') as handle:
            handle.write('#rname\\tstartpos\\tendpos\\tnumreads\\tcovbases\\tcoverage\\tmeandepth\\tmeanbaseq\\tmeanmapq\\n')
            handle.write(f'chromosome\\t1\\t1000\\t10\\t900\\t90\\t{depth}\\t30\\t40\\n')
    PY
    """
}

workflow {
    GENERATE(params.all_rejected ? Channel.of('none') : Channel.of('both', 'denovo_only', 'ref_only', 'none', 'empty_ref'))
    checker = Channel.value(file("${projectDir}/../bin/check_assembly_coverage.py", checkIfExists: true))
    CHECK_DENOVO_COVERAGE(GENERATE.out.denovo, 'denovo', checker)
    CHECK_REF_COVERAGE(GENERATE.out.reference, 'reference-denovo', checker)
    d_reads = GENERATE.out.denovo.map { id, r1, r2, cov -> tuple(id, r1, r2) }
    r_reads = GENERATE.out.reference.map { id, r1, r2, cov -> tuple(id, r1, r2) }
    d_pass = d_reads.join(CHECK_DENOVO_COVERAGE.out.decision)
        .filter { id, r1, r2, passes -> passes == 'true' }
        .map { id, r1, r2, passes -> tuple(id, r1, r2) }
    r_pass = r_reads.join(CHECK_REF_COVERAGE.out.decision)
        .filter { id, r1, r2, passes -> passes == 'true' }
        .map { id, r1, r2, passes -> tuple(id, r1, r2) }
    DENOVO_CANDIDATE(d_pass, 'denovo')
    REF_CANDIDATE(r_pass, 'reference-denovo')
    candidates = DENOVO_CANDIDATE.out.assembly.join(DENOVO_CANDIDATE.out.report)
        .map { id, asm, quast -> tuple(id, 'denovo', asm, quast) }
        .mix(REF_CANDIDATE.out.assembly.join(REF_CANDIDATE.out.report)
            .map { id, asm, quast -> tuple(id, 'reference-denovo', asm, quast) })
        .groupTuple()
        .map { id, methods, assemblies, reports ->
            def d = methods.indexOf('denovo')
            def r = methods.indexOf('reference-denovo')
            def absent = file("${projectDir}/../assets/NO_FILE")
            tuple(id, d >= 0 ? assemblies[d] : absent, r >= 0 ? assemblies[r] : absent,
                  d >= 0 ? reports[d] : absent, r >= 0 ? reports[r] : absent)
        }
    SELECT_BEST_ASSEMBLY(candidates)
    selected = SELECT_BEST_ASSEMBLY.out.assembly.join(SELECT_BEST_ASSEMBLY.out.method.map { method ->
        tuple(method.name.replace('.method.txt', ''), method)
    }).map { id, fasta, method ->
        def expected = id == 'ref_only' ? 'reference-denovo' : 'denovo'
        assert method.text.trim() == expected
        assert fasta.text == ">contig_${expected}\nACGT\n"
        id
    }.toList().map { ids ->
        def expected = params.all_rejected ? [] : ['both', 'denovo_only', 'empty_ref', 'ref_only']
        assert ids.sort() == expected
        "PASS: correct assembly routing and selection for ${ids}"
    }.view()

    ASSEMBLY_COVERAGE_SUMMARY(CHECK_DENOVO_COVERAGE.out.metrics.mix(CHECK_REF_COVERAGE.out.metrics).collect())
    ASSEMBLY_COVERAGE_SUMMARY.out.summary.map { report ->
        assert report.text.contains('Samples with all requested paths skipped: 1')
        assert report.text.contains('none: ASSEMBLY_SKIPPED')
        'PASS: batch coverage summary includes rejected samples'
    }.view()
}
