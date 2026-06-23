module load samtools/1.21
#get unmapped reads into an unmapped.bam
sample=$1
input_bam=$2
thread=$3

#samtools version
unmapped_bam="$sample.unmapped.bam"
unmapped_fq="$sample.unmapped.fq"
samtools view -@ $thread $input_bam -f4 -o $unmapped_bam
#extract unmapped reads into fastq, run mm2 on the unmapped reads to CHM13_T2T
samtools bam2fq $unmapped_bam|gzip -c > $unmapped_fq.gz
less $unmapped_fq.gz|paste - - - -|awk '{print ">"$1"\n"$2}'|gzip -c > $sample.unmapped.fasta.gz

