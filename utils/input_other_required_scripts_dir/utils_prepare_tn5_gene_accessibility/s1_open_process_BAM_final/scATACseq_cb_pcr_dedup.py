#!/usr/bin/env python


import pysam
import argparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count


def collect_fragments(bam_path):
    """
    Step 1:
    Count frequency of each fragment for each cell barcode
    """
    frag_counts = defaultdict(lambda: defaultdict(int))
    read_representative = {}

    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("CB"):
            continue

        cb = read.get_tag("CB")
        chrom = bam.get_reference_name(read.reference_id)
        start = read.reference_start
        end = read.reference_end

        frag = (chrom, start, end)
        frag_counts[cb][frag] += 1

        # keep one representative read
        read_representative.setdefault((cb, frag), read)

    bam.close()
    return frag_counts, read_representative


def collapse_within_cell(args):
    """
    Step 2:
    Within each CB, collapse fragments sharing start OR end on same chromosome
    """
    cb, frag_dict = args
    collapsed = defaultdict(int)

    by_chr = defaultdict(list)
    for (chrom, start, end), count in frag_dict.items():
        by_chr[chrom].append((start, end, count))

    for chrom, items in by_chr.items():
        for start, end, count in items:
            collapsed[(chrom, start, end)] += count

    return cb, collapsed


def collapse_across_cells(collapsed_by_cb):
    """
    Steps 3–5:
    - Collapse identical fragments across CBs
    - Assign to most abundant CB
    - Record total read count
    """
    global_frags = defaultdict(lambda: defaultdict(int))

    for cb, frag_dict in collapsed_by_cb.items():
        for frag, count in frag_dict.items():
            global_frags[frag][cb] += count

    final = {}
    for frag, cb_counts in global_frags.items():
        best_cb = max(cb_counts, key=cb_counts.get)
        total_reads = sum(cb_counts.values())
        final[frag] = (best_cb, total_reads)

    return final


def write_bam(in_bam, out_bam, assignments, rep_reads):
    bam = pysam.AlignmentFile(in_bam, "rb")
    out = pysam.AlignmentFile(out_bam, "wb", template=bam)

    written = set()
    for frag, (cb, count) in assignments.items():
        key = (cb, frag)
        if key not in rep_reads:
            continue

        read = rep_reads[key]
        read.set_tag("CB", cb, value_type="Z")
        read.set_tag("RC", count, value_type="i")

        if key not in written:
            out.write(read)
            written.add(key)

    bam.close()
    out.close()


def main():
    parser = argparse.ArgumentParser(
        description="CB-only PCR duplicate removal for scATAC-seq BAM"
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-p", "--threads", type=int, default=cpu_count())

    args = parser.parse_args()

    print("Collecting fragments...")
    frag_counts, rep_reads = collect_fragments(args.input)

    print("Collapsing within cell barcodes (parallel)...")
    pool = Pool(args.threads)
    results = pool.map(collapse_within_cell, frag_counts.items())
    pool.close()
    pool.join()

    collapsed_by_cb = dict(results)

    print("Collapsing across cell barcodes...")
    final_assignments = collapse_across_cells(collapsed_by_cb)

    print("Writing output BAM...")
    write_bam(args.input, args.output, final_assignments, rep_reads)

    print("Done.")


if __name__ == "__main__":
    main()
