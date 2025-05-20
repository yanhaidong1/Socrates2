#!/usr/bin/env python


import sys
from collections import defaultdict
from natsort import natsorted

def process_file(input_file):
    hash_map = defaultdict(lambda: defaultdict(int))
    its = 0

    if input_file == '-':
        f = sys.stdin
    else:
        f = open(input_file, 'r')

    try:
        for line in f:
            its += 1
            if its % 1_000_000 == 0:
                print(f" - read {its} reads ...", file=sys.stderr)

            cols = line.strip().split('\t')
            cols[0] = 'chr' + cols[0]

            if 'B73V' in cols[0]:
                cols[0] = cols[0].replace('chrB73V4_ctg', 'chrB73V4ctg')

            try:
                site_id = cols[8]
                cell_id = cols[3]
            except IndexError:
                continue  # skip malformed lines

            hash_map[cell_id][site_id] += 1
    finally:
        if f is not sys.stdin:
            f.close()

    cells_sorted = natsorted(hash_map.keys())
    for cell in cells_sorted:
        for site in hash_map[cell]:
            print(f"{site}\t{cell}\t{hash_map[cell][site]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>", file=sys.stderr)
        sys.exit(1)

    input_filename = sys.argv[1]
    process_file(input_filename)