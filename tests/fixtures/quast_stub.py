"""Minimal QUAST stand-in for the final-assembly routing regression test."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('assembly')
parser.add_argument('--min-contig', type=int, required=True)
parser.add_argument('-t')
parser.add_argument('-o', required=True)
args = parser.parse_args()
sequences = Path(args.assembly).read_text().split('>')[1:]
lengths = [len(''.join(record.splitlines()[1:])) for record in sequences]
eligible = [length for length in lengths if length >= args.min_contig]
if not eligible:
    parser.exit(4, 'QUAST stub received an assembly without usable contigs\n')
outdir = Path(args.o)
outdir.mkdir()
(outdir / 'report.tsv').write_text(f'Assembly\ttest\nTotal length\t{sum(eligible)}\n# contigs\t{len(eligible)}\n')
print('QUAST routing stub completed')
