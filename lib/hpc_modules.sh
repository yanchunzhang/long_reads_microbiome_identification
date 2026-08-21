#!/usr/bin/env bash
###############################################################################
# hpc_modules.sh — portable tool loading.
#
# The pipeline was developed on an Lmod cluster (Mount Sinai Minerva) where
# `module load samtools/1.21` is how tools become available. Those calls used to
# be unconditional, which meant every compute job failed immediately under
# `set -euo pipefail` on any system without Lmod — a container, a laptop, a
# conda env, or a SLURM site with different module names.
#
# Source this file and use load_tool instead:
#
#     source "$(dirname "${BASH_SOURCE[0]}")/../../lib/hpc_modules.sh"
#     load_tool samtools samtools/1.21
#     load_tool minimap2 minimap2
#     require_tools samtools minimap2
#
# Resolution order for each tool:
#   1. already on PATH  -> use it, do nothing (container / conda / system)
#   2. Lmod available   -> try `module load`, ignore failure
#   3. otherwise        -> no-op; require_tools reports it clearly
###############################################################################

# Is an environment-module system available in this shell?
_has_lmod() {
    [ -n "${LMOD_CMD:-}" ] && return 0
    type module >/dev/null 2>&1 && return 0
    return 1
}

# load_tool <binary> [module-name]
# Never fails: a missing module must not kill the job before require_tools can
# produce a useful message.
load_tool() {
    local tool="$1" mod="${2:-$1}"
    command -v "$tool" >/dev/null 2>&1 && return 0
    _has_lmod || return 0
    module load "$mod" >/dev/null 2>&1 || true
    return 0
}

# require_tools <binary>...
# Fail fast, naming every missing tool at once, instead of dying mid-pipe with
# "command not found" from inside a multi-stage shell pipeline.
require_tools() {
    local missing=() t
    for t in "$@"; do
        command -v "$t" >/dev/null 2>&1 || missing+=("$t")
    done
    if [ "${#missing[@]}" -gt 0 ]; then
        {
            echo "ERROR: required tool(s) not found on PATH: ${missing[*]}"
            echo
            echo "Provide them in any one of these ways:"
            echo "  - an environment module system (Lmod), e.g. module load samtools"
            echo "  - a conda environment on PATH"
            echo "  - a container built from nextflow/environment.yml"
            echo
            echo "Run check_env.sh to see exactly what is and is not visible to you."
        } >&2
        exit 127
    fi
}
