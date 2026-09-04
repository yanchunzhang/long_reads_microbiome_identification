#!/usr/bin/env bash
###############################################################################
# verify_baseline.sh
#
# Compare a completed test run against the recorded baseline in test/expected/.
#
#   bash test/verify_baseline.sh <run_dir>/CO4_T1
#
# Confirms the pipeline produced the RIGHT answer, not merely that it finished.
# Only order-stable summary outputs are compared: raw BLAST and the chunked
# intermediates can legitimately differ by chunk count and thread count.
#
# Baseline recorded 2026-09-02 from pipeline commit ca6b823, on Minerva, with
# NCBI nt (2024-08-31) and KrakenUniq MicrobialDB (2023-08-08). Different
# database versions WILL produce different results -- that is expected, not a
# regression. See the note at the bottom.
###############################################################################
set -uo pipefail

DIR="${1:-}"
if [[ -z "$DIR" || ! -d "$DIR" ]]; then
    echo "usage: bash test/verify_baseline.sh <run_dir>/CO4_T1" >&2; exit 2
fi
EXPECTED="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/expected"

PASS=0; FAIL=0
ok()  { printf '  \033[32m[ OK ]\033[0m %s\n' "$*"; PASS=$((PASS+1)); }
bad() { printf '  \033[31m[FAIL]\033[0m %s\n' "$*"; FAIL=$((FAIL+1)); }

echo "run      : $DIR"
echo "baseline : $EXPECTED"
echo

echo "Summary outputs"
for f in median_l_adj.txt human_median_length.tsv microbiome.sum_by_length_per_genus.txt \
         blast.ont_adapter_filter.audit.tsv blast.ont_adapter_filtered_out.txt; do
    got="$DIR/CO4_T1.$f"; want="$EXPECTED/CO4_T1.$f"
    if [[ ! -f "$got" ]]; then bad "$f: not produced by the run"
    elif diff -q "$want" "$got" >/dev/null 2>&1; then ok "$f"
    else
        bad "$f: differs from baseline"
        diff "$want" "$got" 2>/dev/null | head -10 | sed 's/^/         /'
    fi
done

echo
echo "Record counts"
while IFS=$'\t' read -r name want; do
    [[ "$name" == "file" ]] && continue
    f="$DIR/$name"
    if [[ ! -f "$f" ]]; then bad "$name: missing"; continue; fi
    case "$name" in
        *.gz)    got=$(zcat "$f" | grep -c '^>') ;;
        *.fasta) got=$(grep -c '^>' "$f") ;;
        *)       got=$(wc -l < "$f") ;;
    esac
    [[ "$got" == "$want" ]] && ok "$name: $got" || bad "$name: got $got, expected $want"
done < "$EXPECTED/record_counts.tsv"

echo
echo "Summary: passed $PASS, failed $FAIL"
if [[ "$FAIL" -gt 0 ]]; then
    cat >&2 <<'EOM'

Differences are not automatically a bug. Check these first:
  * Database versions. NCBI nt grows and its taxonomy changes; the baseline used
    the 2024-08-31 release. A different nt will change assignments legitimately.
  * filter_ont_adapters. The baseline has it enabled (the default). With it
    disabled, blast.microbiome.txt keeps 63 records instead of 62.
If neither explains it, the pipeline has changed behaviour -- investigate.
EOM
    exit 1
fi
echo "Matches the baseline."
