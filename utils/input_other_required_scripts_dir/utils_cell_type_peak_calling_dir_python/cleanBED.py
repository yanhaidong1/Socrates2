#!/usr/bin/env python

import sys


input_file = sys.argv[3]

# Read chromosome lengths from the first input file
chr_lengths = {}


with open(sys.argv[1], 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        chr_lengths[parts[0]] = int(parts[1])

# Compute the adjustment factor from the second argument
factor = float(sys.argv[2]) / 1000000
sys.stderr.write(f" - adjusted factor = {factor}\n")

# Process the second file and apply adjustments


if input_file == '-':
    f = sys.stdin
else:
    f = open(input_file, 'r')


#with open(sys.stdin, 'r') as f:
for line in f:
    parts = line.strip().split('\t')
    chrom = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    value = float(parts[3])

    # Adjust end if it exceeds chromosome length
    if end > chr_lengths.get(chrom, 0):
        end = chr_lengths[chrom] - 1
        if (end - start) < 1:
            continue

    if start >= chr_lengths.get(chrom, 0):
        continue

    adjusted_value = value / factor
    print(f"{chrom}\t{start}\t{end}\t{adjusted_value}")