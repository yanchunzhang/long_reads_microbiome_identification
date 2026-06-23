#!/usr/bin/env bash
###############################################################################
# run_snakemake.sh
#
# Submit the long_reads_microbiome_identification Snakemake workflow to LSF.
#
# This reproduces how the processed datasets were generated, e.g.
#   .../data/processed/datasets/inhouse_mouse/gfm_mouse_mock_spikein_batch1/
#       ONT/long_read_microbiome_fresh_GRCm39/   (run 2026-06-12)
#
# Snakemake is launched with the run directory as the working dir (--directory),
# so it discovers samples from the per-sample subdirs there (each <sample>/
# containing <sample>.bam) and writes outputs back into those subdirs.
# Per-rule cluster resources come from this repo's cluster.yaml.
#
# Usage:
#   bash run_snakemake.sh <run_dir> [options]
#
#   <run_dir>   directory containing per-sample subdirs (<sample>/<sample>.bam).
#               Outputs and logs/ are written here.
#
# Options:
#   -r, --host-ref PATH   minimap2 .mmi for the 2nd host-depletion pass
#                         (overrides config t2t_ref). For MOUSE data this is
#                         GRCm39; for HUMAN data leave unset to use CHM13 T2T.
#   -j, --jobs N          max concurrent cluster jobs        (default: 80)
#   -A, --account ACC     LSF project/account                (default: acc_schzrnas)
#   -n, --dry-run         snakemake -n (print DAG, submit nothing)
#   -u, --unlock          snakemake --unlock then exit
#       --extra "ARGS"    extra args passed verbatim to snakemake
#   -h, --help            show this help
#
# Examples:
#   # Mouse mock run (GRCm39 second pass), as used for batch1:
#   bash run_snakemake.sh \
#     /sc/arion/projects/schzrnas/zhangy40/intratumor_bacteria/data/processed/datasets/\
#inhouse_mouse/gfm_mouse_mock_spikein_batch1/ONT/long_read_microbiome_fresh_GRCm39 \
#     -r /sc/arion/projects/schzrnas/zhangy40/ref/GRCm39/Mus_musculus.GRCm39.dna_sm.toplevel.fa.mmi
#
#   # Dry run first:
#   bash run_snakemake.sh <run_dir> -r <ref.mmi> -n
###############################################################################
set -euo pipefail

# --- environment: put snakemake (and tools) on PATH ------------------------
# The snakefile's shell.prefix only affects RULE jobs on compute nodes; the
# launcher itself (this script) needs snakemake on PATH too.
module load anaconda3 seqkit >/dev/null 2>&1 || true
export PATH="/hpc/users/zhangy40/bin:/hpc/users/zhangy40/schzrnas/softwares/conda/env/myenv/bin:$PATH"

# --- repo locations (this script lives in <repo>/snakemake/) ---------------
REPO_SMK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="$REPO_SMK_DIR/snakefile"
CONFIGFILE="$REPO_SMK_DIR/config.yaml"
CLUSTERYAML="$REPO_SMK_DIR/cluster.yaml"

# --- defaults --------------------------------------------------------------
JOBS=80
ACCOUNT="acc_schzrnas"
HOST_REF=""
DRYRUN=""
UNLOCK=0
EXTRA=""

# --- parse args ------------------------------------------------------------
RUNDIR="${1:-}"
if [[ -z "$RUNDIR" || "$RUNDIR" == "-h" || "$RUNDIR" == "--help" ]]; then
    sed -n '2,55p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 0
fi
shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--host-ref) HOST_REF="$2"; shift 2 ;;
        -j|--jobs)     JOBS="$2"; shift 2 ;;
        -A|--account)  ACCOUNT="$2"; shift 2 ;;
        -n|--dry-run)  DRYRUN="-n"; shift ;;
        -u|--unlock)   UNLOCK=1; shift ;;
        --extra)       EXTRA="$2"; shift 2 ;;
        -h|--help)     sed -n '2,55p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "ERROR: unknown option '$1'" >&2; exit 1 ;;
    esac
done

# --- validate --------------------------------------------------------------
RUNDIR="$(realpath "$RUNDIR")"
[[ -d "$RUNDIR" ]]       || { echo "ERROR: run dir not found: $RUNDIR" >&2; exit 1; }
[[ -f "$SNAKEFILE" ]]    || { echo "ERROR: snakefile not found: $SNAKEFILE" >&2; exit 1; }
[[ -f "$CLUSTERYAML" ]]  || { echo "ERROR: cluster.yaml not found: $CLUSTERYAML" >&2; exit 1; }
if [[ -n "$HOST_REF" ]]; then
    HOST_REF="$(realpath "$HOST_REF")"
    [[ -f "$HOST_REF" ]] || { echo "ERROR: host ref .mmi not found: $HOST_REF" >&2; exit 1; }
fi

# cluster.yaml uses log: "log/{rule}.{wildcards}.%J.log" (relative to workdir),
# and the snakefile creates logs/. Make sure both exist in the run dir.
mkdir -p "$RUNDIR/log" "$RUNDIR/logs"

# --- config overrides (CLI --config wins over config.yaml) -----------------
# scriptsdir MUST be absolute: we run with --directory <run_dir>, so the config's
# default relative "scripts" would resolve inside the run dir (which has none).
CONFIG_OVERRIDES=("scriptsdir=$REPO_SMK_DIR/scripts")
[[ -n "$HOST_REF" ]] && CONFIG_OVERRIDES+=("t2t_ref=$HOST_REF")

# --- LSF submission command (placeholders filled per-rule from cluster.yaml) -
CLUSTER_CMD="bsub -P $ACCOUNT -q {cluster.queue} -n {threads} -W {cluster.time} {cluster.extra} -o {cluster.log} -e {cluster.log}"

echo "=========================================================="
echo " run dir   : $RUNDIR"
echo " snakefile : $SNAKEFILE"
echo " config    : $CONFIGFILE"
echo " cluster   : $CLUSTERYAML  (account=$ACCOUNT)"
echo " host ref  : ${HOST_REF:-<config default: CHM13 T2T>}"
echo " jobs      : $JOBS   dry-run: ${DRYRUN:-no}"
echo "=========================================================="

# --- unlock mode -----------------------------------------------------------
if [[ "$UNLOCK" -eq 1 ]]; then
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIGFILE" \
        --directory "$RUNDIR" --unlock
    echo "unlocked $RUNDIR"; exit 0
fi

# --- launch ----------------------------------------------------------------
# NOTE: run this from a persistent session (tmux/screen) or submit the launcher
# itself as a small long-walltime LSF job, since it stays alive managing the DAG.
snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --directory "$RUNDIR" \
    ${CONFIG_OVERRIDES:+--config "${CONFIG_OVERRIDES[@]}"} \
    --cluster-config "$CLUSTERYAML" \
    --cluster "$CLUSTER_CMD" \
    --jobs "$JOBS" \
    --rerun-incomplete \
    --latency-wait 60 \
    --keep-going \
    --printshellcmds \
    $DRYRUN $EXTRA
