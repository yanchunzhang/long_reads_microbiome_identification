#unmapped_analysis.sh; including steps of get_unmapped reads from input bam file; 2nd round of mapping to t2t ref; get unmapped reads from t2t mapping result bam file.
module load samtools/1.21
module load minimap2

#samtools version
#1st round of getting unmapped reads into an unmapped.bam
sample=$1
input_bam=$2
thread=$3
scriptsdir=$4
t2t_ref=$5

sh $scriptsdir/get_unmapped.sh $sample $input_bam $thread

unmapped_fq="$sample.unmapped.fq.gz"
outprefix="$sample.unmapped.t2t"
sh $scriptsdir/long_read.mm2.no_sort.sh $unmapped_fq $t2t_ref $outprefix $thread
sh $scriptsdir/get_unmapped.sh $sample.after_t2t $sample.unmapped.t2t.mm2.bam $thread


