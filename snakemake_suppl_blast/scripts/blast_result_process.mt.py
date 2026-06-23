#!/usr/bin/env python3

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import islice


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process BLAST hits to find best target per read, with multiprocessing."
    )
    parser.add_argument("--input", required=True, help="Input BLAST result file")
    parser.add_argument("--output", required=True, help="Output summary file")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker processes")
    parser.add_argument(
        "--batch-reads",
        type=int,
        default=5000,
        help="Number of read groups submitted per worker batch"
    )
    return parser.parse_args()


def merge_intervals(intervals):
    if not intervals:
        return 0

    fixed = []
    for s, e in intervals:
        if s <= e:
            fixed.append((s, e))
        else:
            fixed.append((e, s))

    fixed.sort()
    total = 0
    cur_s, cur_e = fixed[0]

    for s, e in fixed[1:]:
        if s <= cur_e + 1:
            if e > cur_e:
                cur_e = e
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e

    total += cur_e - cur_s + 1
    return total


def parse_blast_line(line):
    f = line.rstrip("\n").split("\t")
    return {
        "read_id": f[0],
        "target": f[1],
        "score": float(f[2]),
        "cover1": float(f[3]),
        "cover2": float(f[4]),
        "qstart": int(f[5]),
        "qend": int(f[6]),
        "taxon": f[11],
    }


def choose_best_target(read_hits):
    best_target_name = read_hits[0]["target"]
    best_target_taxon = read_hits[0]["taxon"]

    min_score = min(h["score"] for h in read_hits)
    df_best = [h for h in read_hits if h["score"] == min_score]

    if len({h["target"] for h in df_best}) > 1:
        target_to_intervals = defaultdict(list)
        for h in df_best:
            target_to_intervals[h["target"]].append((h["qstart"], h["qend"]))

        target_lengths = {
            target: merge_intervals(intervals)
            for target, intervals in target_to_intervals.items()
        }

        max_len = max(target_lengths.values())
        best_targets = [t for t, l in target_lengths.items() if l == max_len]

        if len(best_targets) > 1:
            if best_target_name not in best_targets:
                for h in read_hits:
                    if h["target"] in best_targets:
                        best_target_name = h["target"]
                        best_target_taxon = h["taxon"]
                        break
        else:
            if best_targets[0] != best_target_name:
                best_target_name = best_targets[0]
                for h in read_hits:
                    if h["target"] == best_target_name:
                        best_target_taxon = h["taxon"]
                        break

    hit_ranges = [
        (h["qstart"], h["qend"])
        for h in read_hits
        if h["target"] == best_target_name
    ]
    hit_len_combine = merge_intervals(hit_ranges)

    return [read_hits[0]["read_id"], best_target_name, best_target_taxon, hit_len_combine]


def read_groups(blast_file):
    current_read = None
    current_hits = []

    with open(blast_file) as f:
        for line in f:
            if not line.strip():
                continue
            rec = parse_blast_line(line)

            if current_read is None:
                current_read = rec["read_id"]

            if rec["read_id"] != current_read:
                yield current_hits
                current_hits = [rec]
                current_read = rec["read_id"]
            else:
                current_hits.append(rec)

    if current_hits:
        yield current_hits


def chunked(iterable, size):
    it = iter(iterable)
    while True:
        batch = list(islice(it, size))
        if not batch:
            break
        yield batch


def process_batch(batch):
    return [choose_best_target(read_hits) for read_hits in batch]


def main():
    args = parse_arguments()

    with open(args.output, "w") as out:
        if args.threads <= 1:
            for read_hits in read_groups(args.input):
                result = choose_best_target(read_hits)
                out.write("\t".join(map(str, result)) + "\n")
        else:
            with ProcessPoolExecutor(max_workers=args.threads) as ex:
                for results in ex.map(
                    process_batch,
                    chunked(read_groups(args.input), args.batch_reads),
                    chunksize=1
                ):
                    for result in results:
                        out.write("\t".join(map(str, result)) + "\n")


if __name__ == "__main__":
    main()