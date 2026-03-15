#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def human_median_from_stats(stats_file):
    lengths = []
    counts = []

    with open(stats_file) as f:
        for line in f:
            if not line.startswith("RL"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            lengths.append(float(parts[1]))
            counts.append(int(parts[2]))

    df = pd.DataFrame({"length": lengths, "count": counts}).sort_values("length")
    df["cum"] = df["count"].cumsum()
    total = df["count"].sum()

    if total % 2:
        target = (total + 1) // 2
        return df.loc[df["cum"] >= target, "length"].iloc[0]
    else:
        t1 = total // 2
        t2 = t1 + 1
        m1 = df.loc[df["cum"] >= t1, "length"].iloc[0]
        m2 = df.loc[df["cum"] >= t2, "length"].iloc[0]
        return (m1 + m2) / 2


def extract_genus(tax_string):
    if pd.isna(tax_string):
        return None
    for part in str(tax_string).split(";"):
        part = part.strip()
        if part.startswith("g__"):
            return part.replace("g__", "", 1)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", required=True)
    parser.add_argument("--microbe", required=True)
    parser.add_argument("--sum_out", required=True)
    parser.add_argument("--gt5_out", required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    human_median = human_median_from_stats(args.stats)

    if (not os.path.exists(args.microbe)) or os.path.getsize(args.microbe) == 0:
        open(args.sum_out, "w").close()
        open(args.gt5_out, "w").close()
        return

    df = pd.read_csv(
        args.microbe,
        sep="\t",
        header=None,
        usecols=[4, 6]
    )

    if df.empty:
        open(args.sum_out, "w").close()
        open(args.gt5_out, "w").close()
        return

    data = pd.DataFrame({
        "length": pd.to_numeric(df[4], errors="coerce"),
        "genus": df[6].apply(extract_genus)
    }).dropna(subset=["length", "genus"])

    if data.empty:
        open(args.sum_out, "w").close()
        open(args.gt5_out, "w").close()
        return

    summary = (
        data.groupby("genus")["length"]
        .agg(["sum", "count", "median", "mean"])
        .reset_index()
    )
    summary["sample"] = args.sample

    summary.sort_values("sum", ascending=False)[
        ["genus", "sum", "count", "sample"]
    ].to_csv(args.sum_out, sep="\t", header=False, index=False)

    gt5 = summary[summary["count"] >= 5].copy()
    gt5["human_median"] = human_median
    gt5["median_ratio"] = gt5["median"] / gt5["human_median"]

    gt5[
        ["sample", "human_median", "genus", "median", "median_ratio", "mean", "count"]
    ].to_csv(
        args.gt5_out, sep="\t", header=False, index=False
    )


if __name__ == "__main__":
    main()