#!/usr/bin/env python

import sys

hash_dict = {}
its = 0

#input_file = sys.argv[1]

#with open(input_file, 'r') as f:
for line in sys.stdin:
    line = line.strip('\n')
    its += 1
    if its % 1000000 == 0:
        print(f" - read {its} reads ... ", file=sys.stderr)

    line = line.rstrip('\n')
    col = line.split('\t')

    col[0] = 'chr' + col[0]
    if "B73V" in col[0]:
        col[0] = col[0].replace('chrB73V4_ctg', 'chrB73V4ctg')

    id_ = "_".join([col[0], col[6], col[7]])

    if col[3] not in hash_dict:
        hash_dict[col[3]] = {}

    if id_ in hash_dict[col[3]]:
        continue
    else:
        print(f"{id_}\t{col[3]}\t1")
        hash_dict[col[3]][id_] = 1