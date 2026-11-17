nextflow.enable.dsl = 2
params.quast_min_contig = 500
params.outdir = 'final-assembly-test-output'
params.all_rejected = false

include { CHECK_FINAL_ASSEMBLY; FINAL_ASSEMBLY_SUMMARY } from '../modules/local/final_assembly_check'
include { QUAST as QUAST_FINAL } from '../modules/local/quast'

process GENERATE {
    input:
    val sample_id
    output:
    tuple val(sample_id), path("${sample_id}.fasta"), path("${sample_id}.method.txt")
    script:
    def lengths = [boundary: [500], mixed: [600, 30], short: [499, 499], empty: []][sample_id]
    """
    python3 - <<'PY'
    with open('${sample_id}.fasta', 'w') as handle:
        for i, length in enumerate(${lengths}):
            handle.write(f'>contig{i}\\n' + 'A' * length + '\\n')
    with open('${sample_id}.method.txt', 'w') as handle:
        handle.write('denovo\\n')
    PY
    """
}

workflow {
    GENERATE(params.all_rejected ? Channel.of('short', 'empty') : Channel.of('boundary', 'mixed', 'short', 'empty'))
    candidates = GENERATE.out.map { id, fasta, method -> tuple(id, fasta) }
    CHECK_FINAL_ASSEMBLY(candidates, Channel.value(file("${projectDir}/../bin/check_final_assembly.py", checkIfExists: true)))
    passing = candidates.join(CHECK_FINAL_ASSEMBLY.out.decision)
        .filter { id, fasta, passes -> passes == 'true' }
        .map { id, fasta, passes -> tuple(id, fasta) }
    methods = GENERATE.out.map { id, fasta, method -> tuple(id, method) }
        .join(passing.map { id, fasta -> tuple(id) })
        .map { id, method -> id }
    expected = params.all_rejected ? [] : ['boundary', 'mixed']
    methods.toList().map { ids ->
        assert ids.sort() == expected
        'PASS: normal method summaries exclude rejected final assemblies'
    }.view()
    QUAST_FINAL(passing)
    QUAST_FINAL.out.report.map { id, report -> id }.toList().map { ids ->
        assert ids.sort() == expected
        'PASS: only eligible final assemblies reach QUAST'
    }.view()
    FINAL_ASSEMBLY_SUMMARY(CHECK_FINAL_ASSEMBLY.out.metrics.collect())
    FINAL_ASSEMBLY_SUMMARY.out.errors.map { errors ->
        def lines = errors.text.readLines()
        assert lines.size() == 3
        assert lines.any { it.startsWith('short\tSKIPPED\tNO_CONTIGS_ABOVE_MIN_LENGTH') }
        assert lines.any { it.startsWith('empty\tSKIPPED\tNO_CONTIGS_ABOVE_MIN_LENGTH') }
        'PASS: batch errors include short and empty assemblies'
    }.view()
}
