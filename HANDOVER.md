# Handover: running this pipeline without zhangy40's account

This document exists because the pipeline was developed under one account
(`zhangy40`) and originally read several paths out of that account's home
directory, which is mode `0700` and disappears when the account is
deprovisioned. Those paths have been repointed; this is what remains to know.

**Outside Mount Sinai entirely?** This document is about the Sinai handover.
See the "Running outside Mount Sinai" section of [README.md](README.md) for
scheduler profiles, the config template, and tool installation.

---

## 1. The two things you need access to

| What | Where | Size | Fallback if you lack access |
|---|---|---|---|
| **Tools** (conda env, taxonkit, taxdump, this repo) | `/sc/arion/projects/schzrnas/zhangy40/softwares/` | ~3.4 GB | **Yes** — copy staged at `/sc/arion/projects/fangg03a/long_read_microbiome/` |
| **Databases** (BLAST `nt`, 2 KrakenUniq DBs, CHM13 T2T) | `/sc/arion/projects/schzrnas/zhangy40/` | **~1.2 TB** | **No** — you must be in group `schzrnas` |

The databases are the binding constraint. They are far too large to duplicate,
so **anyone running this pipeline needs membership in the `schzrnas` Unix
group.** Request it from Minerva HPC support; the `schzrnas` allocation is owned
by `fangg03` (the PI), so this is an in-lab approval.

Note that the Unix group (`schzrnas`) is a *separate thing* from the LSF billing
account (`acc_schzrnas`). You may well have one and not the other — see §4.

## 2. Quick start

```bash
# 1. Preflight — verifies every DB, reference and tool is readable BY YOU.
#    Takes seconds. Do this before submitting anything.
bash check_env.sh

# 2. Stage your data as <run_dir>/<sample>/<sample>.bam
mkdir -p /sc/arion/projects/<your_proj>/$USER/run1/SAMPLE1
cp SAMPLE1.bam /sc/arion/projects/<your_proj>/$USER/run1/SAMPLE1/

# 3. Dry run, then real run
bash snakemake/run_snakemake.sh /sc/arion/projects/<your_proj>/$USER/run1 -n
bash snakemake/run_snakemake.sh /sc/arion/projects/<your_proj>/$USER/run1
```

Outputs are written back into each `<sample>/` subdirectory. Nothing is written
into `zhangy40`'s space, so each person runs against their own run directory.

For **mouse** data add `-r <GRCm39.mmi>`; human data uses the CHM13 T2T default.

## 3. Environment

There is nothing to install. Tool locations resolve automatically at run time:

1. `$PIPELINE_ENV_BIN` / `$PIPELINE_TOOLS` if set, otherwise
2. the canonical `schzrnas` paths if readable, otherwise
3. the `fangg03a` copy.

To force the copy explicitly:

```bash
export PIPELINE_ENV_BIN=/sc/arion/projects/fangg03a/long_read_microbiome/env/myenv/bin
export PIPELINE_TOOLS=/sc/arion/projects/fangg03a/long_read_microbiome/bin
```

**Do not run `conda activate`.** The pipeline only prepends the env's `bin` to
`PATH`. This matters: `module load anaconda3` pulls in miniforge3, which
prepends its own `bin` and will shadow the env's python 3.7 with python 3.13
(losing pandas 1.3.5). Modules must be loaded *before* the `PATH` export —
`snakefile`, `run_snakemake.sh` and `check_env.sh` all do this in the right
order, so just don't reorder them.

The env supplies only **krakenuniq 1.0.4, taxonkit 0.14.1, python 3.7.12,
pandas 1.3.5, snakemake 7.24.0**. samtools, minimap2, blastn and seqkit come
from `module load`.

### Rebuilding the environment

Prefer using the shared env as-is — it is the environment that produced the
published results. If you must rebuild:

```bash
conda create -p ./myenv --file env/myenv.explicit.lock.txt   # exact, 341 packages
```

`nextflow/environment.yml` is a human-readable equivalent, but solving from it
may drift (python 3.7 is EOL). **Do not let pandas float to 2.x** — its groupby
and median semantics differ, and the MLA metric depends on both.

`nextflow/pipeline.def` has never been built successfully: its `%post` uses
`apt-get`, which needs fakeroot, and Minerva users are not in `/etc/subuid`. It
also pins no versions. The `-profile singularity` profile is therefore untested;
build and validate against a known sample before trusting it.

## 4. LSF billing account

`run_snakemake.sh` defaults to `acc_schzrnas`. If you bill to a different
account your jobs are **rejected at submission**, which is confusing because
`check_env.sh` will have passed. Either:

```bash
bash snakemake/run_snakemake.sh <run_dir> -A acc_fangg03a
# or
export LSF_ACCOUNT=acc_fangg03a
```

Nextflow: `--lsf_project acc_fangg03a`, or the same `LSF_ACCOUNT` variable.

## 5. Long-running launcher

The Snakemake launcher stays alive managing the DAG for the whole run. Start it
inside `tmux`/`screen`, or submit it as its own long-walltime LSF job. If it
dies mid-run, `run_snakemake.sh <run_dir> -u` unlocks the directory.

## 6. Paths that were changed for the handover

For anyone auditing against older runs — these were repointed, and all resolve
to identical content:

| Setting | Was | Now |
|---|---|---|
| `blastdb` | `/hpc/users/zhangy40/schzrnas/softwares/blast_db/nt` | `/sc/arion/projects/schzrnas/zhangy40/softwares/blast_db/nt` |
| `taxonkit_data_dir` | `/hpc/users/zhangy40/.taxonkit` | `.../softwares/taxdump` (byte-identical, verified with `cmp`) |
| `scratch_dir` | `/sc/arion/scratch/zhangy40` | `/sc/arion/scratch/$USER` |
| `PATH` prefix | `/hpc/users/zhangy40/bin:...` | resolved via `PIPELINE_ENV_BIN` / `PIPELINE_TOOLS` |

The first two were the same files all along — `/hpc/users/zhangy40/schzrnas` is
a symlink into project space. Only the *route* was private.
