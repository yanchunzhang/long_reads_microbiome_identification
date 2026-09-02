#!/usr/bin/env python3
"""Remove ONT barcode/adapter-derived bases from BLAST query coverage.

The input microbiome table is headerless and has the pipeline's seven columns:
read ID, representative subject, taxid, covered bp, read length, query coverage,
and lineage.  A read is not discarded merely because an ONT sequence is found.
Instead, bases where a representative-target BLAST HSP overlaps an ONT hit are
subtracted from covered bp.  The assignment is removed only when the adjusted
query coverage no longer exceeds the configured threshold.
"""

import argparse
import csv
import shutil
from collections import defaultdict


AUDIT_HEADER = [
    "read_id", "representative_subject", "read_length", "ont_elements",
    "ont_read_intervals", "original_covered_bp", "raw_hsp_union_bp",
    "ont_overlap_with_blast_bp", "adjusted_covered_bp", "original_query_coverage",
    "adjusted_query_coverage", "decision", "reason",
]


def merge_intervals(intervals):
    """Merge 1-based, closed intervals."""
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def interval_bp(intervals):
    return sum(end - start + 1 for start, end in intervals)


def overlap_bp(left, right):
    total = 0
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo <= hi:
            total += hi - lo + 1
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--microbiome", required=True)
    parser.add_argument("--raw-blast", required=True)
    parser.add_argument("--ont-hits", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit-output", required=True)
    parser.add_argument("--filtered-output", required=True)
    parser.add_argument("--end-window", type=int, default=150)
    parser.add_argument("--min-overlap", type=int, default=18)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-query-coverage", type=float, default=0.5)
    parser.add_argument("--disabled", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()

    rows = []
    by_read = {}
    with open(args.microbiome) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                raise ValueError("Expected at least 7 tab-separated microbiome columns")
            record = {
                "fields": fields,
                "line": line,
                "subject": fields[1],
                "covered": int(float(fields[3])),
                "length": int(float(fields[4])),
                "qcov": float(fields[5]),
            }
            rows.append((fields[0], record))
            by_read[fields[0]] = record

    if args.disabled:
        shutil.copyfile(args.microbiome, args.output)
        open(args.filtered_output, "w").close()
        with open(args.audit_output, "w", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerow(AUDIT_HEADER)
        return

    # Raw megaBLAST: qseqid, sseqid, evalue, pident, length, qstart, qend, ...
    hsp_intervals = defaultdict(list)
    with open(args.raw_blast) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or fields[0] not in by_read:
                continue
            if fields[1] != by_read[fields[0]]["subject"]:
                continue
            start, end = sorted((int(fields[5]), int(fields[6])))
            hsp_intervals[fields[0]].append((start, end))

    # Adapter BLAST format:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    ont_intervals = defaultdict(list)
    ont_elements = defaultdict(set)
    with open(args.ont_hits) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[0] not in by_read:
                continue
            read_id, element = fields[0], fields[1]
            identity, aligned = float(fields[2]), int(fields[3])
            start, end = sorted((int(fields[6]), int(fields[7])))
            read_len = by_read[read_id]["length"]
            near_end = start <= args.end_window or end > read_len - args.end_window
            if aligned >= args.min_overlap and identity >= args.min_identity and near_end:
                ont_intervals[read_id].append((start, end))
                ont_elements[read_id].add(element)

    with open(args.output, "w") as retained, \
         open(args.filtered_output, "w") as filtered, \
         open(args.audit_output, "w", newline="") as audit_handle:
        audit = csv.writer(audit_handle, delimiter="\t", lineterminator="\n")
        audit.writerow(AUDIT_HEADER)

        for read_id, record in rows:
            if read_id not in ont_intervals:
                retained.write(record["line"])
                continue

            adapter_union = merge_intervals(ont_intervals[read_id])
            hsp_union = merge_intervals(hsp_intervals.get(read_id, []))
            overlap = overlap_bp(hsp_union, adapter_union)
            adjusted_bp = max(0, record["covered"] - overlap)
            adjusted_qcov = adjusted_bp / record["length"] if record["length"] else 0.0
            keep = adjusted_qcov > args.min_query_coverage

            fields = list(record["fields"])
            fields[3] = str(adjusted_bp)
            fields[5] = "{:.10g}".format(adjusted_qcov)
            adjusted_line = "\t".join(fields) + "\n"
            if keep:
                retained.write(adjusted_line)
                decision = "retain"
                reason = ("ONT sequence detected but adjusted BLAST query coverage passes"
                          if overlap else "ONT sequence does not overlap representative BLAST support")
            else:
                filtered.write(adjusted_line)
                decision = "filter"
                reason = "adjusted BLAST query coverage does not exceed threshold"

            audit.writerow([
                read_id, record["subject"], record["length"],
                ";".join(sorted(ont_elements[read_id])),
                ";".join("{}-{}".format(start, end) for start, end in adapter_union),
                record["covered"], interval_bp(hsp_union), overlap, adjusted_bp,
                "{:.10g}".format(record["qcov"]), "{:.10g}".format(adjusted_qcov),
                decision, reason,
            ])


if __name__ == "__main__":
    main()
