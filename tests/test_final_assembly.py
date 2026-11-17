import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('final_assembly', Path(__file__).resolve().parents[1] / 'bin/check_final_assembly.py')
final_check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(final_check)


class FinalAssemblyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.assembly = Path(self.temp.name) / 'assembly.fasta'

    def assess(self, lengths, minimum=500):
        with self.assembly.open('w') as handle:
            for i, length in enumerate(lengths):
                handle.write(f'>contig{i}\n' + 'A' * length + '\n')
        return final_check.assess('sample', self.assembly, minimum)

    def test_boundary_is_inclusive(self):
        self.assertEqual(self.assess([499])['status'], 'SKIPPED')
        self.assertEqual(self.assess([500])['status'], 'PASS')

    def test_combined_short_contigs_do_not_pass(self):
        result = self.assess([499, 499])
        self.assertEqual(result['reason'], 'NO_CONTIGS_ABOVE_MIN_LENGTH')
        self.assertEqual(result['total_bases'], 998)
        self.assertEqual(result['longest_contig'], 499)

    def test_one_eligible_contig_is_enough(self):
        result = self.assess([20, 600, 100])
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['eligible_contigs'], 1)
        self.assertEqual(result['eligible_bases'], 600)

    def test_empty_and_header_only_assemblies_are_skipped(self):
        self.assertEqual(self.assess([])['status'], 'SKIPPED')
        self.assertEqual(self.assess([0])['status'], 'SKIPPED')

    def test_multiline_sequence_lengths_are_combined(self):
        self.assembly.write_text('>contig\n' + 'A' * 250 + '\n' + 'C' * 250 + '\n')
        self.assertEqual(final_check.assess('sample', self.assembly, 500)['longest_contig'], 500)

    def test_sequence_without_header_is_an_input_error(self):
        self.assembly.write_text('ACGT\n')
        with self.assertRaises(ValueError):
            final_check.assess('sample', self.assembly, 500)


if __name__ == '__main__':
    unittest.main()
