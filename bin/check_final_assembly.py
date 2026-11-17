#!/usr/bin/env python3
"""Reject final assemblies without any contigs meeting QUAST's length cutoff."""

import argparse
import csv
from pathlib import Path


def assess(sample, assembly, minimum_length):
    lengths = []
    current = None
    with open(assembly) as handle:
        for line in handle:
            if line.startswith('>'):
                if current is not None:
                    lengths.append(current)
                current = 0
            elif line.strip():
                if current is None:
                    raise ValueError(f'Sequence before the first FASTA header in {assembly}')
                current += len(''.join(line.split()))
    if current is not None:
        lengths.append(current)
    eligible = [length for length in lengths if length >= minimum_length]
    return {
        'sample': sample,
        'status': 'PASS' if eligible else 'SKIPPED',
        'reason': 'PASS' if eligible else 'NO_CONTIGS_ABOVE_MIN_LENGTH',
        'minimum_contig_length': minimum_length,
        'total_contigs': len(lengths),
        'total_bases': sum(lengths),
        'longest_contig': max(lengths, default=0),
        'eligible_contigs': len(eligible),
        'eligible_bases': sum(eligible),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--assembly', required=True)
    parser.add_argument('--minimum-length', type=int, required=True)
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--error-report', required=True)
    args = parser.parse_args()
    if args.minimum_length <= 0:
        parser.error('--minimum-length must be positive')
    result = assess(args.sample, args.assembly, args.minimum_length)
    with open(args.metrics, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result), delimiter='\t')
        writer.writeheader()
        writer.writerow(result)
    if result['status'] != 'PASS':
        Path(args.error_report).write_text(
            'CtGAP Final Assembly Error Report\n'
            '=================================\n'
            + ''.join(f'{key}: {value}\n' for key, value in result.items())
            + '\nNo contigs meet the minimum length required by QUAST.\n'
            'Final QUAST and downstream analyses were skipped for this sample.\n'
            'Inspect the filtering report and the unfiltered assembly.\n'
        )
    print('true' if result['status'] == 'PASS' else 'false')


if __name__ == '__main__':
    main()
