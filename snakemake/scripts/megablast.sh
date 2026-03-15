ml blast/2.13.0+
query=$1
db=$2
out=$3
thread=$4

blastn -query $query -task megablast -db $db -out $out -outfmt "6 qseqid sseqid evalue pident length qstart qend sstart send bitscore stitle staxids sscinames" -max_target_seqs 5 -num_threads $thread
