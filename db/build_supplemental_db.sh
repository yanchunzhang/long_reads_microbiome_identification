#!/bin/bash
#BSUB -P acc_schzrnas
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=64000] span[hosts=1]"
#BSUB -W 72:00
#BSUB -J build_suppl_krakenuniq_db
#BSUB -o logs/build_supplemental_db_%J.out
#BSUB -e logs/build_supplemental_db_%J.err

# Build a supplemental KrakenUniq database to fill coverage gaps in the primary
# kuniq_microbialdb_minus_kdb.20230808 database.
#
# What this adds:
#   - RefSeq viral genomes (Complete Genome only)
#     Primary DB has most viruses but gaps include Influenza A/B, HCV, HBV, HSV-2
#   - High-quality RefSeq fungal genomes (Complete/Chromosome level assemblies)
#     beyond the EuPathDB subset already in the primary DB
#
# Quality criteria:
#   Viruses  : assembly_level = Complete Genome only (Salzberg Sci Transl Med 2025)
#              no min contig length — complete viral genomes are full-length by definition
#   Fungi    : assembly_level ∈ {Complete Genome, Chromosome} (PathSeq-T2T Cell 2026)
#              genome_rep = Full
#              minimum contig length = 5,000 bp (PathSeq standard for non-viruses)
#
# Usage:
#   bsub < db/build_supplemental_db.sh
#
# After building, run KrakenUniq twice — once per DB — and merge outputs.
# See db/supplemental_database.md for pipeline integration details.

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
SUPDB="/sc/arion/projects/schzrnas/zhangy40/softwares/kuniq_supplemental_vf_db"
THREADS=8
MIN_CONTIG_FUNGI=5000       # bp; PathSeq standard for non-viral contigs
WORK_DIR="/sc/arion/scratch/zhangy40/suppl_db_build"
LOG_PREFIX="[$(date '+%Y-%m-%d %H:%M:%S')]"

# ── Environment ───────────────────────────────────────────────────────────────
module load anaconda3
set +u  # conda activate scripts reference unbound vars (libxml2 known issue)
conda activate myenv
set -u

mkdir -p "$SUPDB" "$WORK_DIR"

echo "$LOG_PREFIX Starting supplemental KrakenUniq database build"
echo "$LOG_PREFIX Output DB  : $SUPDB"
echo "$LOG_PREFIX Work dir   : $WORK_DIR"
echo "$LOG_PREFIX Threads    : $THREADS"

# ── Step 1: Download NCBI taxonomy ────────────────────────────────────────────
# Must be done first; provides the taxid → lineage mappings used by all
# subsequent krakenuniq-build steps.
echo "$LOG_PREFIX [Step 1/4] Downloading NCBI taxonomy..."
if [[ ! -f "$SUPDB/taxDB" ]]; then
    krakenuniq-download --db "$SUPDB" taxonomy
else
    echo "$LOG_PREFIX   taxDB already present, skipping."
fi

# ── Shared helper: FASTA contig filter ───────────────────────────────────────
# Reformats headers to |kraken:taxid|TAXID and drops contigs below min_len.
# Written once; reused for both viral and fungal downloads.
FILTER_PY="$WORK_DIR/filter_contigs.py"
cat > "$FILTER_PY" << 'PYEOF'
#!/usr/bin/env python3
"""Filter FASTA by minimum contig length and add kraken:taxid header tags."""
import sys
import gzip
import re

in_gz   = sys.argv[1]
out_fa  = sys.argv[2]
taxid   = sys.argv[3]
min_len = int(sys.argv[4])

opener = gzip.open if in_gz.endswith('.gz') else open

kept = dropped = 0
with opener(in_gz, 'rt') as fh, open(out_fa, 'w') as out:
    header, seq_parts = None, []
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if header is not None:
                seq = ''.join(seq_parts)
                if len(seq) >= min_len:
                    seqid = re.sub(r'\s.*', '', header[1:])
                    out.write(f'>{seqid}|kraken:taxid|{taxid}\n{seq}\n')
                    kept += 1
                else:
                    dropped += 1
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        seq = ''.join(seq_parts)
        if len(seq) >= min_len:
            seqid = re.sub(r'\s.*', '', header[1:])
            out.write(f'>{seqid}|kraken:taxid|{taxid}\n{seq}\n')
            kept += 1
        else:
            dropped += 1

print(f'  kept={kept} dropped_short={dropped}', file=sys.stderr)
PYEOF
chmod +x "$FILTER_PY"

# ── Shared download-filter-add loop ──────────────────────────────────────────
# Used for both viral and fungal sections below.
download_and_add() {
    local summary_file="$1"
    local min_len="$2"
    local label="$3"
    local TOTAL=0 SKIPPED=0 FAILED=0

    while IFS=$'\t' read -r -a fields; do
        ftp_path="${fields[19]}"
        taxid="${fields[5]}"
        asm_name=$(basename "$ftp_path")
        done_flag="$WORK_DIR/${asm_name}.done"

        if [[ -f "$done_flag" ]]; then
            (( SKIPPED++ )) || true
            continue
        fi

        fasta_url="${ftp_path}/${asm_name}_genomic.fna.gz"
        local_gz="$WORK_DIR/${asm_name}_genomic.fna.gz"
        local_filt="$WORK_DIR/${asm_name}.filtered.fna"

        if ! wget -q -O "$local_gz" "$fasta_url"; then
            echo "WARN: download failed for $fasta_url"
            (( FAILED++ )) || true
            rm -f "$local_gz"
            continue
        fi

        python3 "$FILTER_PY" "$local_gz" "$local_filt" "$taxid" "$min_len"

        # krakenuniq-build --add-to-library exits 255 even on success (Perl quirk in v1.0.4).
        # Use || true so set -e does not abort the loop; success is confirmed by the library file.
        krakenuniq-build --add-to-library "$local_filt" --db "$SUPDB" || true

        touch "$done_flag"
        rm -f "$local_gz" "$local_filt"
        (( TOTAL++ )) || true

    done < "$summary_file"

    echo "$LOG_PREFIX   $label: added=$TOTAL  skipped(already done)=$SKIPPED  failed=$FAILED"
}

# ── Step 2: RefSeq viral genomes (Complete Genome only) ──────────────────────
# Salzberg (Sci Transl Med 2025) uses "all RefSeq viral complete genomes" —
# complete only, not chromosome-level. This is the stricter of the two paper
# standards and appropriate for viruses since genuine complete viral genomes
# are inherently full-length single sequences.
# No minimum contig length — complete viral genomes need no further filtering.
# Fills gaps in the primary DB: Influenza A/B, HCV, HBV, HSV-2, etc.
echo "$LOG_PREFIX [Step 2/4] Downloading RefSeq viral genomes (Complete Genome only)..."

VIRAL_SUMMARY="$WORK_DIR/assembly_summary_viral.txt"
VIRAL_FILTERED="$WORK_DIR/assembly_summary_viral.filtered.txt"

if [[ ! -f "$VIRAL_FILTERED" ]]; then
    wget -q -O "$VIRAL_SUMMARY" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt"

    awk -F'\t' \
        'NR>2 && $12=="Complete Genome" && $20!="na"' \
        "$VIRAL_SUMMARY" > "$VIRAL_FILTERED"

    N=$(wc -l < "$VIRAL_FILTERED")
    echo "$LOG_PREFIX   Retained $N viral assemblies after quality filter."
fi

download_and_add "$VIRAL_FILTERED" 1 "Viruses"

# ── Step 3: RefSeq fungal genomes (Complete Genome + Chromosome) ──────────────
# Same assembly-level filter as viruses. Short contigs (<5 kbp) additionally
# dropped — standard for non-viral sequences (PathSeq; Walker et al. 2018).
echo "$LOG_PREFIX [Step 3/4] Downloading RefSeq fungal genomes (Complete + Chromosome)..."

FUNGI_SUMMARY="$WORK_DIR/assembly_summary_fungi.txt"
FUNGI_FILTERED="$WORK_DIR/assembly_summary_fungi.filtered.txt"

if [[ ! -f "$FUNGI_FILTERED" ]]; then
    wget -q -O "$FUNGI_SUMMARY" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt"

    # genome_rep = Full additionally required for fungi (more assemblies available;
    # partial representations add ambiguous k-mers without coverage benefit)
    awk -F'\t' \
        'NR>2 && ($12=="Complete Genome" || $12=="Chromosome") && $14=="Full" && $20!="na"' \
        "$FUNGI_SUMMARY" > "$FUNGI_FILTERED"

    N=$(wc -l < "$FUNGI_FILTERED")
    echo "$LOG_PREFIX   Retained $N fungal assemblies after quality filter."
fi

download_and_add "$FUNGI_FILTERED" "$MIN_CONTIG_FUNGI" "Fungi"

# ── Step 4: Build the database ────────────────────────────────────────────────
# --kmer-len 31 and --minimizer-len 15 match the primary database defaults.
# This step is the most memory- and time-intensive (~64 GB RAM, several hours).
echo "$LOG_PREFIX [Step 4/4] Building KrakenUniq database..."
export JELLYFISH_BIN="${JELLYFISH_BIN:-jellyfish}"  # build_db.sh references this with set -u
krakenuniq-build \
    --build \
    --db "$SUPDB" \
    --threads "$THREADS" \
    --kmer-len 31 \
    --minimizer-len 15

echo "$LOG_PREFIX Build complete."
echo "$LOG_PREFIX Supplemental database path: $SUPDB"
echo "$LOG_PREFIX"
echo "$LOG_PREFIX Next step: update nextflow/nextflow.config to add:"
echo "$LOG_PREFIX   params.kraken_db_suppl = \"$SUPDB\""
