#!/usr/bin/env bash
# unmapped_analysis.sh
# Get unmapped reads from an input BAM, map them to the T2T reference, then
# emit reads still unmapped after T2T filtering.

set -euo pipefail

module load samtools/1.21
module load minimap2

sample=$1
input_bam=$(realpath "$2")
outdir=$(realpath "$3")
thread=$4
scriptsdir=$(realpath "$5")
t2t_ref=$(realpath "$6")

mkdir -p "$outdir"

(
    cd "$outdir"

    sh "$scriptsdir/get_unmapped.sh" "$sample" "$input_bam" "$thread"

    unmapped_fq="$sample.unmapped.fq.gz"
    outprefix="$sample.unmapped.t2t"
    sh "$scriptsdir/long_read.mm2.no_sort.sh" "$unmapped_fq" "$t2t_ref" "$outprefix" "$thread"

    sh "$scriptsdir/get_unmapped.sh" "$sample.after_t2t" "$sample.unmapped.t2t.mm2.bam" "$thread"
)
