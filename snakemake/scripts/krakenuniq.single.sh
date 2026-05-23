#!/usr/bin/env bash
# krakenuniq.single.sh
# Run KrakenUniq on a single FASTA input.
# Environment is expected to be set up by the caller (shell.prefix / beforeScript).
#
# Usage: krakenuniq.single.sh <DB> <input_fasta> <prefix> <threads>

set -euo pipefail

DB=$1
input=$2
prefix=$3
thread=$4

classified_out="${prefix}.classified_by_krakenuniq"

krakenuniq --preload --db "$DB" --threads "$thread" \
    --classified-out "$classified_out" \
    --report "${prefix}.krakenuniq.report" \
    "$input" > "${prefix}.krakenuniq"
