ml samtools
input=$1
thread=$2
samtools view -@ $thread $input -h -F4| samtools stats > $input.mapped_human_reads_only.stats.txt
