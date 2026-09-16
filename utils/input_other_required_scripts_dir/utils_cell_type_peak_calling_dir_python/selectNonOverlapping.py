#!/usr/bin/env python

import sys

def parse_coord(entry):
    """Split a string like chr_100_200_5 into parts."""
    chrom, start, end, score = entry.split('_')
    return chrom, int(start), int(end), float(score)

def is_overlap(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)

input_file = sys.argv[1]

#with open(input_file, 'r') as f:
#    for line in sys.stdin:
#        line = line.strip('\n')
for line in sys.stdin:
    line = line.strip('\n')
    if "B73" in line:
        line = line.replace("B73V4_ctg", "B73V4ctg")

    cols = line.split('\t')
    chrom = cols[0]
    start = int(cols[1])
    end = int(cols[2])
    field = cols[3]

    if ',' not in field:
        # Single entry
        parts = field.split('_')
        score = parts[-1]
        print(f"{chrom}\t{start}\t{end}\t{score}")
    else:
        # Multiple comma-separated entries
        sites = field.split(',')
        site_dict = {}
        for site in sites:
            _, s, e, score = parse_coord(site)
            site_dict[site] = float(score)

        # Sort by score descending
        sorted_sites = sorted(site_dict.keys(), key=lambda x: -site_dict[x])

        taken = []
        for i, site in enumerate(sorted_sites):
            _, s1, e1, score = parse_coord(site)
            overlap = False
            for taken_site in taken:
                _, ts1, te1, _ = parse_coord(taken_site)
                if is_overlap(s1, e1, ts1, te1):
                    overlap = True
                    break
            if not overlap:
                print(f"{chrom}\t{s1}\t{e1}\t{score}")
                taken.append(site)