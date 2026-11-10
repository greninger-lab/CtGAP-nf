import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('assembly_coverage', Path(__file__).resolve().parents[1] / 'bin/check_assembly_coverage.py')
coverage_check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coverage_check)


class AssemblyCoverageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.r1 = self.fastq('r1', 1)
        self.r2 = self.fastq('r2', 1)

    def fastq(self, name, count):
        path = self.root / f'{name}.fastq.gz'
        with gzip.open(path, 'wt') as handle:
            for i in range(count):
                handle.write(f'@read{i}\nACGT\n+\nIIII\n')
        return path

    def coverage(self, depths, sample_column=False):
        path = self.root / 'coverage.tsv'
        prefix = 'sample_name\t' if sample_column else ''
        header = '#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n'
        with path.open('w') as handle:
            handle.write(prefix + header)
            for i, depth in enumerate(depths):
                handle.write(('sample\t' if sample_column else '') + f'ref{i}\t1\t1000\t10\t900\t90\t{depth}\t30\t40\n')
        return path

    def assess(self, depths, minimum=10, sample_column=False):
        return coverage_check.assess('sample', 'denovo', self.r1, self.r2, self.coverage(depths, sample_column), minimum)

    def test_depth_threshold_boundary(self):
        self.assertEqual(self.assess([9.99])['reason'], 'LOW_MAPPED_DEPTH')
        self.assertEqual(self.assess([10])['status'], 'PASS')
        self.assertEqual(self.assess([5], minimum=5)['status'], 'PASS')

    def test_reference_set_is_not_diluted_by_unused_references(self):
        result = self.assess([0] * 23 + [12], sample_column=True)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['reference'], 'ref23')
        self.assertEqual(result['mean_depth'], 12)

    def test_empty_extracted_reads_override_high_alignment_depth(self):
        self.r1 = self.fastq('empty1', 0)
        self.r2 = self.fastq('empty2', 0)
        self.assertEqual(self.assess([50])['reason'], 'NO_READ_PAIRS')
        self.assertEqual(self.assess([0], minimum=0)['status'], 'SKIPPED')

    def test_unbalanced_mates_are_skipped(self):
        self.r2 = self.fastq('extra', 2)
        self.assertEqual(self.assess([50])['reason'], 'UNBALANCED_READ_PAIRS')

    def test_header_only_coverage_is_zero(self):
        self.assertEqual(self.assess([])['mean_depth'], 0)
        self.assertEqual(self.assess([])['reason'], 'LOW_MAPPED_DEPTH')

    def test_invalid_coverage_is_not_silently_treated_as_low_depth(self):
        for value in ['NaN', '-1']:
            with self.assertRaises(ValueError):
                self.assess([value])

    def test_corrupt_fastq_is_not_silently_treated_as_empty(self):
        with gzip.open(self.r1, 'wt') as handle:
            handle.write('@read\nACGT\n+\nI\n')
        with self.assertRaises(ValueError):
            self.assess([50])


if __name__ == '__main__':
    unittest.main()
