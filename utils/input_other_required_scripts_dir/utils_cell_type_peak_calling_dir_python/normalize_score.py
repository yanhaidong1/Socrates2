#!/usr/bin/env python

import sys

input_file = sys.argv[1]

# First pass: compute the total of column 5 (index 4)
total = 0
with open(input_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        total += float(parts[4])

# Compute the scale factor
scale_factor = total / 1000000

# Second pass: scale and output
with open(input_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        score = float(parts[4]) / scale_factor
        print(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{score}")
