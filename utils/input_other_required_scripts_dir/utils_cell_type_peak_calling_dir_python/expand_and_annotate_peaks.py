#!/usr/bin/env python

import sys

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    cols = line.split("\t")
    chrom = cols[0]
    start = int(cols[1]) - 250
    end = int(cols[2]) + 250
    score = cols[3]

    if start >= 0:
        name = f"{chrom}_{start}_{end}_{score}"
        print(f"{chrom}\t{start}\t{end}\t{score}\t{name}")