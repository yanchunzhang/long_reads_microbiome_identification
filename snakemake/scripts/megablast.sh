#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/hpc_modules.sh"
load_tool blastn blast/2.13.0+
require_tools blastn

query=$1
db=$2
out=$3
thread=$4

blastn -query $query -task megablast -db $db -out $out -outfmt "6 qseqid sseqid evalue pident length qstart qend sstart send bitscore stitle staxids sscinames" -max_target_seqs 5 -max_hsps 10 -num_threads $thread
