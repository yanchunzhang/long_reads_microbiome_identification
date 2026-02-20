import pandas as pd
import pybedtools
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process BLAST hits to find best targets per read")
    parser.add_argument('--input', required=True, help='Input BLAST result file')
    parser.add_argument('--output', required=True, help='Output summary file')
    parser.add_argument('--n_start', type=int, default=0, help='Start index of reads to process')
    parser.add_argument('--n_to_process', type=int, default=-1, help='Number of reads to process (-1 means all)')
    return parser.parse_args()

def process_read_hits(df, read_id):
    df_flt = df[df[0] == read_id].reset_index(drop=True)
    best_target_name = df_flt.at[0, 1]
    best_target_taxon = df_flt.at[0, 11]

    # Get rows with minimum alignment score (column 2)
    df_best = df_flt[df_flt[2] == df_flt[2].min()]

    if df_best[1].nunique() > 1:
        # Calculate adjusted coverage and merge intervals
        df_best = df_best.assign(cover_adjusted=df_best[3] * df_best[4] / 100)
        df2_bed = df_best[[1, 5, 6]].rename(columns={1: 'chrom', 5: 'start', 6: 'end'})
        merged = pybedtools.BedTool.from_dataframe(df2_bed).sort().merge().to_dataframe(names=['chrom', 'start', 'end'])
        merged['len'] = merged['end'] - merged['start'] + 1

        max_len = merged.groupby('chrom')['len'].sum()
        max_len = max_len[max_len == max_len.max()]

        if len(max_len) > 1:
            if best_target_name not in max_len.index:
                best_target_name = df_flt[df_flt[1].isin(max_len.index)].iloc[0, 1]
                best_target_taxon = df_flt[df_flt[1].isin(max_len.index)].iloc[0, 11]
        elif max_len.index[0] != best_target_name:
            best_target_name = max_len.index[0]
            best_target_taxon = df_flt[df_flt[1] == best_target_name].iloc[0, 11]

    # Combine hit ranges for best target
    hit_ranges = df_flt[df_flt[1] == best_target_name][[1, 5, 6]].rename(columns={1: 'chrom', 5: 'start', 6: 'end'})
    merged_ranges = pybedtools.BedTool.from_dataframe(hit_ranges).sort().merge().to_dataframe(names=['chrom', 'start', 'end'])
    hit_len_combine = (merged_ranges['end'] - merged_ranges['start']).sum() + len(merged_ranges)

    return [read_id, best_target_name, best_target_taxon, hit_len_combine]

def main():
    args = parse_arguments()
    df = pd.read_table(args.input, header=None, sep="\t")
    read_ids = df[0].unique()

    n_end = None if args.n_to_process == -1 else args.n_start + args.n_to_process
    read_ids_to_process = read_ids[args.n_start:n_end]

    df = df[df[0].isin(read_ids_to_process)]
    #df = pd.DataFrame()
    with open(args.output, "w") as f:
        for read_id in read_ids_to_process:
            result = process_read_hits(df, read_id)
            f.write("\t".join(str(x) for x in result) + "\n")

if __name__ == "__main__":
    main()
