#!/usr/bin/env bash

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/lib/hpc_modules.sh"
load_tool samtools samtools
require_tools samtools

input=$1
thread=$2
samtools view -@ $thread $input -h -F4| samtools stats > $input.mapped_human_reads_only.stats.txt
