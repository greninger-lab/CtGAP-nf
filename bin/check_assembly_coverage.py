#!/usr/bin/env python3
"""Check mapped depth and the actual paired FASTQs supplied to an assembler."""

import argparse
import csv
import gzip
import math
from pathlib import Path


def fastq_stats(path):
    reads = bases = 0
    with gzip.open(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip('\r\n')
            separator = handle.readline()
            quality = handle.readline().rstrip('\r\n')
            if not header.startswith('@') or not separator.startswith('+') or not sequence or len(sequence) != len(quality):
                raise ValueError(f'Invalid FASTQ record in {path} at read {reads + 1}')
            reads += 1
            bases += len(sequence)
    return reads, bases


def best_coverage(path):
    best = ('none', 0.0, 0.0)
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        reference_key = '#rname' if '#rname' in fields else 'rname'
        if reference_key not in fields or 'meandepth' not in fields or 'coverage' not in fields:
            raise ValueError(f'Missing samtools coverage columns in {path}')
        for row in reader:
            depth = float(row['meandepth'])
            breadth = float(row['coverage'])
            if not math.isfinite(depth) or depth < 0 or not math.isfinite(breadth) or not 0 <= breadth <= 100:
                raise ValueError(f'Invalid coverage metrics in {path}')
            if best[0] == 'none' or depth > best[1]:
                best = (row[reference_key], depth, breadth)
    return best


def assess(sample, method, r1, r2, coverage, minimum_depth):
    n1, bases1 = fastq_stats(r1)
    n2, bases2 = fastq_stats(r2)
    reference, depth, breadth = best_coverage(coverage)
    if n1 == 0 or n2 == 0:
        reason = 'NO_READ_PAIRS'
    elif n1 != n2:
        reason = 'UNBALANCED_READ_PAIRS'
    elif depth < minimum_depth:
        reason = 'LOW_MAPPED_DEPTH'
    else:
        reason = 'PASS'
    return {
        'sample': sample, 'assembly_path': method, 'status': 'PASS' if reason == 'PASS' else 'SKIPPED',
        'reason': reason, 'reference': reference, 'mean_depth': depth, 'coverage_percent': breadth,
        'minimum_depth': minimum_depth, 'r1_reads': n1, 'r2_reads': n2, 'read_bases': bases1 + bases2,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--method', required=True, choices=['denovo', 'reference-denovo'])
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', required=True)
    parser.add_argument('--coverage', required=True)
    parser.add_argument('--minimum-depth', required=True, type=float)
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--error-report', required=True)
    args = parser.parse_args()
    if not math.isfinite(args.minimum_depth) or args.minimum_depth < 0:
        parser.error('--minimum-depth must be finite and non-negative')
    result = assess(args.sample, args.method, args.r1, args.r2, args.coverage, args.minimum_depth)
    with open(args.metrics, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result), delimiter='\t')
        writer.writeheader()
        writer.writerow(result)
    if result['status'] != 'PASS':
        Path(args.error_report).write_text(
            'CtGAP Assembly Coverage Error Report\n'
            '====================================\n'
            + ''.join(f'{key}: {value}\n' for key, value in result.items())
            + '\nThis assembly path was skipped before Shovill.\n'
            'In auto mode, the other path may still produce an assembly.\n'
        )
    print('true' if result['status'] == 'PASS' else 'false')


if __name__ == '__main__':
    main()
