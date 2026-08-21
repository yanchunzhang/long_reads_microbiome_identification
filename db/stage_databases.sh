#!/usr/bin/env bash
###############################################################################
# stage_databases.sh
#
# Copy the reference databases into a second location so they can be read by a
# group that lacks access to the originals.
#
# ONLY RUN THIS IF THE GROUP-MEMBERSHIP REQUEST HAS FAILED. Adding the readers
# to the `schzrnas` Unix group costs one ticket and zero bytes; this costs
# ~1.2 TB and creates a second copy that will drift from the original. Prefer
# the ticket.
#
# This does NOT move anything. The originals stay where they are: their paths
# are referenced by ~236 files across the intratumor_bacteria project (146 for
# CHM13 alone) and by other members of the schzrnas allocation.
#
# Usage:
#   bash stage_databases.sh [dest_dir] [--dry-run] [--only NAME]
#
#   dest_dir   default /sc/arion/projects/fangg03a/zhangy40/long_read_microbiome/db
#              This is the location the pipeline already probes as a fallback,
#              so a completed copy is picked up with no config change.
#
#   --only     stage just one: kraken | suppl | blast | chm13 | taxdump
#
# Safe to re-run: rsync resumes rather than restarting.
###############################################################################
set -euo pipefail

SRC=/sc/arion/projects/schzrnas/zhangy40
DEST="${1:-/sc/arion/projects/fangg03a/zhangy40/long_read_microbiome/db}"
[[ "${1:-}" == --* ]] && DEST=/sc/arion/projects/fangg03a/zhangy40/long_read_microbiome/db
DRYRUN=""; ONLY=""
for a in "$@"; do
    case "$a" in
        --dry-run) DRYRUN="--dry-run" ;;
        --only)    ONLY="NEXT" ;;
        *)         [[ "$ONLY" == "NEXT" ]] && ONLY="$a" ;;
    esac
done

# name|source|approx size
ITEMS=(
  "taxdump|$SRC/softwares/taxdump|465M"
  "kraken|$SRC/softwares/kuniq_microbialdb_minus_kdb.20230808|544G"
  "suppl|$SRC/softwares/kuniq_supplemental_vf_db_v2|130G"
  "blast|$SRC/softwares/blast_db|534G"
  "chm13|$SRC/ref/CHM13_T2T|28G"
)

echo "=========================================================="
echo " source : $SRC"
echo " dest   : $DEST"
echo " mode   : ${DRYRUN:-copy}${ONLY:+   only=$ONLY}"
echo "=========================================================="
echo
echo "Free space at destination:"; df -h "$(dirname "$DEST")" | tail -1
echo
echo "NOTE: df shows FILESYSTEM free space, not your allocation's quota."
echo "      Confirm the fangg03a quota with HPC before copying ~1.2 TB."
echo

if [[ -z "$DRYRUN" ]]; then
    read -r -p "Copy ~1.2 TB now? Type 'yes' to continue: " ans
    [[ "$ans" == "yes" ]] || { echo "aborted"; exit 1; }
fi

mkdir -p "$DEST"
for item in "${ITEMS[@]}"; do
    IFS='|' read -r name src size <<< "$item"
    [[ -n "$ONLY" && "$ONLY" != "$name" ]] && continue
    tgt="$DEST/$(basename "$src")"
    echo "── $name ($size) → $tgt"
    if [[ ! -r "$src" ]]; then echo "   SKIP: source not readable"; continue; fi
    rsync -a --info=progress2 $DRYRUN "$src/" "$tgt/"
done

if [[ -z "$DRYRUN" ]]; then
    # CRITICAL: rsync -a preserves the SOURCE group. Because the operator is in
    # both groups it can set schzrnas on the copy, which would leave the whole
    # point defeated -- the files unreadable by exactly the people they are for.
    # This bit the tools copy: 88,758 files landed with the wrong group.
    echo "── fixing group ownership and permissions"
    chgrp -R fangg03a "$DEST"
    chmod -R g+rX "$DEST"
    find "$DEST" -type d -exec chmod g+s {} \; 2>/dev/null || true

    echo
    echo "── verification"
    bad_group=$(find "$DEST" ! -group fangg03a 2>/dev/null | wc -l)
    unreadable=$(find "$DEST" ! -perm -g+r 2>/dev/null | wc -l)
    printf '   wrong group: %s   not group-readable: %s\n' "$bad_group" "$unreadable"
    [[ "$bad_group" -eq 0 && "$unreadable" -eq 0 ]] \
        && echo "   OK" || { echo "   *** FIX BEFORE USE ***"; exit 1; }
    echo
    echo "Done. The pipeline probes this location automatically -- run"
    echo "  bash check_env.sh"
    echo "as a user WITHOUT schzrnas access to confirm the fallback resolves."
fi
