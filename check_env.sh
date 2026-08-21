#!/usr/bin/env bash
###############################################################################
# check_env.sh
#
# Preflight for the long_reads_microbiome_identification pipeline.
#
# Verifies that every database, reference, and tool the pipeline needs is
# actually READABLE BY YOU before any job is submitted. Run this first on a new
# account: a missing permission surfaces here in seconds instead of as a failed
# LSF job an hour into the queue.
#
# Usage:
#   bash check_env.sh [run_dir]
#
#   run_dir   optional; if given, also checks that you can write there and that
#             it contains per-sample subdirs of the form <sample>/<sample>.bam
#
# Exit status: 0 = all good, 1 = at least one blocking problem.
###############################################################################
set -uo pipefail

PASS=0; FAIL=0; WARN=0
ok()   { printf '  \033[32m[ OK ]\033[0m %s\n' "$*"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31m[FAIL]\033[0m %s\n' "$*"; FAIL=$((FAIL+1)); }
warn() { printf '  \033[33m[WARN]\033[0m %s\n' "$*"; WARN=$((WARN+1)); }
hdr()  { printf '\n\033[1m%s\033[0m\n' "$*"; }

# --- resolve site paths exactly as run_snakemake.sh does ------------------
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

SCHZ=/sc/arion/projects/schzrnas/zhangy40
FANG=/sc/arion/projects/fangg03a/zhangy40/long_read_microbiome

ENV_BIN="$(resolve_site_path PIPELINE_ENV_BIN "$SCHZ/softwares/conda/env/myenv/bin" "$FANG/env/myenv/bin")"
TOOLS_DIR="$(resolve_site_path PIPELINE_TOOLS "$SCHZ/softwares" "$FANG/bin")"

hdr "Identity"
echo "  user   : $(id -un)  (uid $(id -u))"
echo "  groups : $(id -Gn)"
if id -Gn | tr ' ' '\n' | grep -qx schzrnas; then
    ok "in group 'schzrnas' — canonical paths available"
    IN_SCHZ=1
else
    warn "NOT in group 'schzrnas' — falling back to the fangg03a copy for tools."
    warn "     The 1.2 TB databases live in schzrnas space and have no fallback;"
    warn "     if the DB checks below fail, request schzrnas group membership."
    IN_SCHZ=0
fi

hdr "Resolved tool locations"
echo "  PIPELINE_ENV_BIN : $ENV_BIN"
echo "  PIPELINE_TOOLS   : $TOOLS_DIR"
[[ -r "$ENV_BIN"   && -x "$ENV_BIN"   ]] && ok "env bin readable"   || bad "env bin NOT readable: $ENV_BIN"
[[ -r "$TOOLS_DIR" && -x "$TOOLS_DIR" ]] && ok "tools dir readable" || bad "tools dir NOT readable: $TOOLS_DIR"

# ORDER MATTERS: load modules FIRST, then prepend our PATH. `module load
# anaconda3` pulls in miniforge3, which prepends its own bin and would otherwise
# shadow myenv's python 3.7 with python 3.13 (and lose pandas 1.3.5).
# This mirrors the ordering in snakefile shell.prefix and run_snakemake.sh.
module load anaconda3 seqkit >/dev/null 2>&1 || warn "'module load anaconda3 seqkit' failed"
module load blast/2.13.0+     >/dev/null 2>&1 || true
module load samtools/1.21     >/dev/null 2>&1 || true
module load minimap2/2.24     >/dev/null 2>&1 || true
export PATH="$ENV_BIN:$TOOLS_DIR:$PATH"

hdr "Tools on PATH"
check_tool() {  # name, version-cmd, expected-substring (optional)
    local name="$1" vcmd="$2" want="${3:-}"
    local path; path="$(command -v "$name" 2>/dev/null)"
    if [[ -z "$path" ]]; then bad "$name: NOT FOUND"; return; fi
    local ver; ver="$(eval "$vcmd" 2>&1 | head -1 | tr -d '\r')"
    if [[ -n "$want" && "$ver" != *"$want"* ]]; then
        warn "$name: $ver  (expected to contain '$want')  [$path]"
    else
        ok "$name: ${ver:-present}  [$path]"
    fi
}
check_tool krakenuniq "krakenuniq --version" "1.0.4"
check_tool taxonkit   "taxonkit version"     "0.14.1"
check_tool python     "python --version"     "3.7"
check_tool samtools   "samtools --version"
check_tool minimap2   "minimap2 --version"
check_tool blastn     "blastn -version"
check_tool seqkit     "seqkit version"
check_tool snakemake  "snakemake --version"

if command -v python >/dev/null 2>&1; then
    pv="$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null)"
    if [[ -z "$pv" ]]; then bad "pandas: not importable"
    elif [[ "$pv" == 1.3.5 ]]; then ok "pandas: $pv"
    else warn "pandas: $pv — published results used 1.3.5; 2.x changes groupby/median behaviour"
    fi
fi

hdr "Databases and references"
check_path() {  # label, path, [one expected file inside]
    local label="$1" p="$2" probe="${3:-}"
    if [[ ! -e "$p" ]]; then bad "$label: does not exist ($p)"; return; fi
    if [[ ! -r "$p" ]]; then bad "$label: EXISTS BUT NOT READABLE ($p)"; return; fi
    if [[ -n "$probe" && ! -r "$p/$probe" ]]; then bad "$label: missing/unreadable $probe in $p"; return; fi
    ok "$label ($(du -sh "$p" 2>/dev/null | cut -f1))"
}
check_path "kraken MicrobialDB"  "$SCHZ/softwares/kuniq_microbialdb_minus_kdb.20230808" "taxDB"
check_path "kraken suppl VF DB"  "$SCHZ/softwares/kuniq_supplemental_vf_db_v2"          "taxDB"
check_path "CHM13 T2T index"     "$SCHZ/ref/CHM13_T2T/chm13v2.0.mmi"

TAXDUMP="$(resolve_site_path PIPELINE_TAXONKIT_DB "$SCHZ/softwares/taxdump" "$FANG/db/taxdump")"
for f in nodes.dmp names.dmp merged.dmp delnodes.dmp; do
    [[ -r "$TAXDUMP/$f" ]] && ok "taxdump/$f" || bad "taxdump/$f unreadable (dir: $TAXDUMP)"
done

BLASTDB_DIR="$SCHZ/softwares/blast_db"
if [[ -r "$BLASTDB_DIR" ]] && compgen -G "$BLASTDB_DIR/nt.*" >/dev/null; then
    ok "BLAST nt ($(ls "$BLASTDB_DIR"/nt.*.nin 2>/dev/null | wc -l) volumes)"
    if command -v blastdbcmd >/dev/null 2>&1; then
        if blastdbcmd -db "$BLASTDB_DIR/nt" -info >/dev/null 2>&1; then ok "BLAST nt opens cleanly"
        else bad "BLAST nt present but blastdbcmd cannot open it"; fi
    fi
else
    bad "BLAST nt unreadable: $BLASTDB_DIR"
fi

hdr "LSF"
if command -v bsub >/dev/null 2>&1; then
    ok "bsub present"
    ACC="${LSF_ACCOUNT:-}"
    if [[ -n "$ACC" ]]; then echo "  LSF_ACCOUNT=$ACC"
    elif [[ "$IN_SCHZ" == 1 ]]; then echo "  will default to acc_schzrnas"
    else warn "not in schzrnas: set LSF_ACCOUNT (e.g. acc_fangg03a) or pass -A, or jobs will be rejected"
    fi
else
    warn "bsub not found — are you on a login node?"
fi

hdr "Scratch"
SCRATCH="/sc/arion/scratch/$(id -un)"
if [[ -d "$SCRATCH" && -w "$SCRATCH" ]]; then ok "scratch writable: $SCRATCH"
else warn "scratch missing or not writable: $SCRATCH (megablast/blast_process stage there)"; fi

if [[ $# -ge 1 ]]; then
    hdr "Run directory: $1"
    if [[ ! -d "$1" ]]; then bad "not a directory"
    elif [[ ! -w "$1" ]]; then bad "not writable by you"
    else
        ok "exists and is writable"
        n=$(find "$1" -mindepth 2 -maxdepth 2 -name '*.bam' 2>/dev/null | wc -l)
        [[ "$n" -gt 0 ]] && ok "found $n sample BAM(s) as <sample>/<sample>.bam" \
                         || warn "no <sample>/<sample>.bam found — nothing to run"
    fi
fi

hdr "Summary"
printf '  passed %d   warnings %d   failures %d\n' "$PASS" "$WARN" "$FAIL"
if [[ "$FAIL" -gt 0 ]]; then
    printf '\n\033[31mNOT ready to run.\033[0m Fix the [FAIL] items above.\n'
    exit 1
fi
printf '\n\033[32mReady to run.\033[0m Next: bash snakemake/run_snakemake.sh <run_dir> -n\n'
