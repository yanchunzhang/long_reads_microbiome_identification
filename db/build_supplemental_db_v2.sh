#!/bin/bash
#BSUB -P acc_schzrnas
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=64000] span[hosts=1]"
#BSUB -W 72:00
#BSUB -J build_suppl_krakenuniq_db
#BSUB -o logs/build_supplemental_db_%J.out
#BSUB -e logs/build_supplemental_db_%J.err

# Supplemental KrakenUniq database — clean build using krakenuniq-download.
#
# Replaces build_supplemental_db.sh which used a custom NCBI FTP download loop
# and krakenuniq-build --add-to-library. That approach had a Perl bug in v1.0.4
# that silently skipped seqid2taxid.map creation. krakenuniq-download handles
# download, header tagging, .map file creation, and contig filtering natively.
#
# What this adds:
#   - RefSeq viral genomes (Complete Genome only)
#   - RefSeq fungal genomes (Complete Genome + Chromosome, genome_rep=Full)
#
# Quality criteria match published standards:
#   Viruses : Complete Genome only (Salzberg Sci Transl Med 2025)
#   Fungi   : Complete Genome or Chromosome, genome_rep=Full (PathSeq-T2T Cell 2026)
#             minimum contig length = 5,000 bp (PathSeq standard)
#
# Usage:
#   bsub < db/build_supplemental_db_v2.sh
#
# NOTE: Run on an empty SUPDB directory. Downloads are resumable — existing
# files are skipped unless --force is passed to krakenuniq-download.

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
SUPDB="/sc/arion/projects/schzrnas/zhangy40/softwares/kuniq_supplemental_vf_db_v2"
THREADS=8
MIN_CONTIG_FUNGI=5000

# ── Environment ───────────────────────────────────────────────────────────────
module load anaconda3
set +u
conda activate myenv
set -u

export JELLYFISH_BIN="${JELLYFISH_BIN:-jellyfish}"

mkdir -p "$SUPDB"
LOG_PREFIX="[$(date '+%Y-%m-%d %H:%M:%S')]"

echo "$LOG_PREFIX Starting supplemental KrakenUniq database build (v2)"
echo "$LOG_PREFIX Output DB : $SUPDB"
echo "$LOG_PREFIX Threads   : $THREADS"

# ── Step 1: Taxonomy ──────────────────────────────────────────────────────────
echo "$LOG_PREFIX [Step 1/4] Downloading NCBI taxonomy..."
if [[ ! -f "$SUPDB/taxDB" ]]; then
    krakenuniq-download --db "$SUPDB" taxonomy
else
    echo "$LOG_PREFIX   taxDB already present, skipping."
fi

# ── Step 2: Viral genomes ─────────────────────────────────────────────────────
# Complete Genome only — Salzberg Sci Transl Med 2025 standard.
# No minimum contig length; complete viral genomes are single full-length sequences.
echo "$LOG_PREFIX [Step 2/4] Downloading RefSeq viral genomes (Complete Genome only)..."
krakenuniq-download \
    --db "$SUPDB" \
    --threads "$THREADS" \
    refseq/viral/Complete_Genome

# ── Step 3: Fungal genomes ────────────────────────────────────────────────────
# Complete Genome + Chromosome, genome_rep=Full — PathSeq-T2T Cell 2026 standard.
# Contig filter (≥5,000 bp) applied by --min-seq-len.
echo "$LOG_PREFIX [Step 3/4] Downloading RefSeq fungal genomes (Complete + Chromosome, genome_rep=Full)..."
krakenuniq-download \
    --db "$SUPDB" \
    --threads "$THREADS" \
    --min-seq-len "$MIN_CONTIG_FUNGI" \
    'refseq/fungi/Complete_Genome/genome_rep=Full'

krakenuniq-download \
    --db "$SUPDB" \
    --threads "$THREADS" \
    --min-seq-len "$MIN_CONTIG_FUNGI" \
    'refseq/fungi/Chromosome/genome_rep=Full'

# ── Step 4: Build ─────────────────────────────────────────────────────────────
echo "$LOG_PREFIX [Step 4/4] Building KrakenUniq database..."
krakenuniq-build \
    --build \
    --db "$SUPDB" \
    --threads "$THREADS" \
    --kmer-len 31 \
    --minimizer-len 15

echo "$LOG_PREFIX Build complete."
echo "$LOG_PREFIX Supplemental database path: $SUPDB"
