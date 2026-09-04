#!/usr/bin/env python3
"""Filter reads dominated by merged ONT barcode/adapter construct sequence.

The input microbiome table is headerless and has the pipeline's seven columns:
read ID, representative subject, taxid, covered bp, read length, query coverage,
and lineage. Qualifying ONT-hit query intervals are merged without double
counting. A read is filtered when their union occupies at least the configured
fraction of the complete read. The ordinary upstream BLAST criteria remain
unchanged; BLAST coverage is not adjusted a second time here.
"""

import argparse
import csv
import shutil
from collections import defaultdict


AUDIT_HEADER = [
    "read_id", "representative_subject", "read_length", "ont_elements",
    "ont_read_intervals", "ont_technical_covered_bp", "ont_technical_fraction",
    "original_covered_bp", "original_query_coverage", "decision", "reason",
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--microbiome", required=True)
    parser.add_argument("--ont-hits", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit-output", required=True)
    parser.add_argument("--filtered-output", required=True)
    parser.add_argument("--min-overlap", type=int, default=12)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-technical-fraction", type=float, default=0.40)
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
            if aligned >= args.min_overlap and identity >= args.min_identity:
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
            technical_bp = interval_bp(adapter_union)
            technical_fraction = (technical_bp / record["length"]
                                  if record["length"] else 0.0)
            should_filter = technical_fraction >= args.min_technical_fraction

            if not should_filter:
                retained.write(record["line"])
                decision = "retain"
                reason = "merged ONT technical coverage is below threshold"
            else:
                filtered.write(record["line"])
                decision = "filter"
                reason = "merged ONT technical coverage meets or exceeds threshold"

            audit.writerow([
                read_id, record["subject"], record["length"],
                ";".join(sorted(ont_elements[read_id])),
                ";".join("{}-{}".format(start, end) for start, end in adapter_union),
                technical_bp, "{:.10g}".format(technical_fraction),
                record["covered"], "{:.10g}".format(record["qcov"]), decision, reason,
            ])


if __name__ == "__main__":
    main()
