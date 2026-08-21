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
#   -A, --account ACC     LSF project/account   (default: $LSF_ACCOUNT or acc_schzrnas)
#                         Members of another lab MUST set this, e.g. acc_fangg03a,
#                         or LSF rejects the job at submission.
#   -l, --local           run on this machine instead of submitting to a
#                         scheduler (use --cores N via --extra). Needed off-LSF.
#       --cluster-cmd C   submit command template, overriding the LSF default.
#                         Also settable as $SNAKEMAKE_CLUSTER_CMD. SLURM e.g.:
#                           'sbatch -p {cluster.queue} -c {threads} -t {cluster.time}'
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

# Print the banner comment block (between the two ### rules), stripped of "# ".
usage() {
    awk 'NR>2 && /^#{10,}/ {exit} NR>2 {sub(/^# ?/, ""); print}' "${BASH_SOURCE[0]}"
}

# --- environment: put snakemake (and tools) on PATH ------------------------
# The snakefile's shell.prefix only affects RULE jobs on compute nodes; the
# launcher itself (this script) needs snakemake on PATH too.
module load anaconda3 seqkit >/dev/null 2>&1 || true

# Resolve tool locations: env var wins, else first readable candidate.
# schzrnas members get the canonical paths; everyone else falls through to the
# fangg03a copy automatically.
resolve_site_path() {
    local var_name="$1"; shift
    local override="${!var_name:-}"
    if [[ -n "$override" ]]; then echo "$override"; return; fi
    local c
    for c in "$@"; do
        if [[ -r "$c" && -x "$c" ]]; then echo "$c"; return; fi
    done
    echo "$1"
}

ENV_BIN="$(resolve_site_path PIPELINE_ENV_BIN \
    /sc/arion/projects/schzrnas/zhangy40/softwares/conda/env/myenv/bin \
    /sc/arion/projects/fangg03a/zhangy40/long_read_microbiome/env/myenv/bin)"
TOOLS_DIR="$(resolve_site_path PIPELINE_TOOLS \
    /sc/arion/projects/schzrnas/zhangy40/softwares \
    /sc/arion/projects/fangg03a/zhangy40/long_read_microbiome/bin)"
export PATH="$ENV_BIN:$TOOLS_DIR:$PATH"

# --- repo locations (this script lives in <repo>/snakemake/) ---------------
REPO_SMK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="$REPO_SMK_DIR/snakefile"
CONFIGFILE="$REPO_SMK_DIR/config.yaml"
CLUSTERYAML="$REPO_SMK_DIR/cluster.yaml"

# --- defaults --------------------------------------------------------------
JOBS=80
ACCOUNT="${LSF_ACCOUNT:-acc_schzrnas}"
HOST_REF=""
DRYRUN=""
UNLOCK=0
EXTRA=""
LOCAL=0
CLUSTER_CMD_OVERRIDE="${SNAKEMAKE_CLUSTER_CMD:-}"

# --- parse args ------------------------------------------------------------
RUNDIR="${1:-}"
if [[ -z "$RUNDIR" || "$RUNDIR" == "-h" || "$RUNDIR" == "--help" ]]; then
    usage
    exit 0
fi
shift
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--host-ref) HOST_REF="$2"; shift 2 ;;
        -j|--jobs)     JOBS="$2"; shift 2 ;;
        -A|--account)  ACCOUNT="$2"; shift 2 ;;
        -l|--local)    LOCAL=1; shift ;;
        --cluster-cmd) CLUSTER_CMD_OVERRIDE="$2"; shift 2 ;;
        -n|--dry-run)  DRYRUN="-n"; shift ;;
        -u|--unlock)   UNLOCK=1; shift ;;
        --extra)       EXTRA="$2"; shift 2 ;;
        -h|--help)     usage; exit 0 ;;
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
# Default submit command is LSF (Minerva). Override with --cluster-cmd or
# $SNAKEMAKE_CLUSTER_CMD for another scheduler, or use --local for no scheduler.
# NOTE: written as an if/else, NOT ${VAR:-default}. The template contains {...}
# placeholders and parameter expansion ends at the first `}`, which silently
# mangles the command into `-q {cluster.queue -n {threads} ...}`.
if [[ -n "$CLUSTER_CMD_OVERRIDE" ]]; then
    CLUSTER_CMD="$CLUSTER_CMD_OVERRIDE"
else
    CLUSTER_CMD="bsub -P $ACCOUNT -q {cluster.queue} -n {threads} -W {cluster.time} {cluster.extra} -o {cluster.log} -e {cluster.log}"
fi

echo "=========================================================="
echo " run dir   : $RUNDIR"
echo " snakefile : $SNAKEFILE"
echo " config    : $CONFIGFILE"
if [[ "$LOCAL" -eq 1 ]]; then
    echo " execution : LOCAL (no scheduler; pass --extra \"--cores N\")"
else
    echo " cluster   : $CLUSTERYAML  (account=$ACCOUNT)"
    echo " submit    : $CLUSTER_CMD"
fi
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
SCHED_ARGS=(--cluster-config "$CLUSTERYAML" --cluster "$CLUSTER_CMD" --jobs "$JOBS")
if [[ "$LOCAL" -eq 1 ]]; then
    # No scheduler: let snakemake run rules on this machine. Pass --cores via
    # --extra to control parallelism, e.g. --extra "--cores 8".
    SCHED_ARGS=()
fi

snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --directory "$RUNDIR" \
    ${CONFIG_OVERRIDES:+--config "${CONFIG_OVERRIDES[@]}"} \
    ${SCHED_ARGS[@]+"${SCHED_ARGS[@]}"} \
    --rerun-incomplete \
    --latency-wait 60 \
    --keep-going \
    --printshellcmds \
    $DRYRUN $EXTRA
