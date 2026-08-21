#!/usr/bin/env bash
#unmapped_analysis.sh; including steps of get_unmapped reads from input bam file; 2nd round of mapping to t2t ref; get unmapped reads from t2t mapping result bam file.
#samtools version
#1st round of getting unmapped reads into an unmapped.bam

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/hpc_modules.sh"
load_tool samtools samtools/1.21
load_tool minimap2 minimap2
require_tools samtools minimap2

sample=$1
input_bam=$2
thread=$3

scriptsdir="scripts"

sh $scriptsdir/get_unmapped.sh $sample $input_bam $thread

t2t_ref="/sc/arion/projects/schzrnas/zhangy40/ref/CHM13_T2T/chm13v2.0.mmi"
unmapped_fq="$sample.unmapped.fq.gz"
outprefix="$sample.unmapped.t2t"
sh $scriptsdir/long_read.mm2.no_sort.sh $unmapped_fq $t2t_ref $outprefix $thread
sh $scriptsdir/get_unmapped.sh $sample.after_t2t $sample.unmapped.t2t.mm2.bam $thread


