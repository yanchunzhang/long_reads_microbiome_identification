#long_read.mm2.no_sort.sh
module load samtools/1.21
module load minimap2

fq=$1
ref=$2
outprefix=$3
thread=$4
type=$5

mm2_ax="map-ont"

if [ "$type" = "pacbio" ];
then
  mm2_ax="map-hifi"
fi

minimap2 -ax $mm2_ax -t $thread $ref $fq > $outprefix.mm2.sam
samtools view $outprefix.mm2.sam -O bam -o $outprefix.mm2.bam && rm $outprefix.mm2.sam

