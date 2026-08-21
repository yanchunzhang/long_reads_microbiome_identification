# long_reads_microbiome_identification

A pipeline for detecting genuine microbiome signals in low-biomass tissues using long-read sequencing (ONT/PacBio). Available in both **Nextflow** and **Snakemake** implementations.

> **Running this on a Mount Sinai account that isn't `zhangy40`?** Start with
> [HANDOVER.md](HANDOVER.md), then run `bash check_env.sh` to verify you can
> actually read every database and tool before submitting any jobs.

## Background

Microbial DNA detection in human tissues is often confounded by contamination. This pipeline implements a DNA fragment length-based metric to discriminate genuine microbial DNA (long fragments) from contaminant DNA (short fragments), as described in:

> Zhang et al. (2026). *Critical assessment of intratumor and low-biomass microbiome using long-read sequencing.* bioRxiv. https://www.biorxiv.org/content/10.64898/2026.02.02.703393v2

## Pipeline Overview

```
BAM (long-read, aligned to human reference)
  │
  ├── 1. Extract unmapped reads (two-round filtering: GRCh38 → T2T CHM13v2)
  │
  ├── 2a. KrakenUniq — primary microbial DB
  ├── 2b. KrakenUniq — supplemental viral+fungal DB (optional, use_suppl_db=true)
  │         └── merge primary + supplemental hits (primary DB has priority)
  │
  ├── 3. megaBLAST validation (parallelized by chunk)
  │
  ├── 4. Fragment length normalization (host vs. microbial read lengths)
  │
  └── 5. Per-sample genus-level microbiome abundance
```

## How It Works

1. **Human read filtering** — Unmapped reads are extracted from the input BAM (hg38-aligned) using samtools, then re-aligned to the T2T CHM13v2 reference with minimap2. Reads remaining unmapped after both steps are retained as candidate non-human reads.
2. **Human read length distribution** — `samtools stats` is run on the original BAM to obtain the median human read length, used later for normalization.
3. **Taxonomic classification** — Candidate reads are classified with KrakenUniq against the primary MicrobialDB database. If `use_suppl_db=true`, a supplemental viral+fungal database is run in parallel and the results are merged (primary DB has priority; duplicate read IDs are deduplicated).
4. **BLAST validation** — KrakenUniq-classified microbial reads are aligned to NCBI nt using megaBLAST. The best hit per read (largest alignment coverage) is retained; only hits with coverage >0.5 are kept.
5. **Genus-level aggregation** — Microbial reads are grouped by genus (≥5 reads per genus required).
6. **Median(L)adj calculation** — The key metric:

   **Median(L)adj** = median microbial read length (per genus) / median human read length

   Genuine microbial DNA has long fragments; contaminant DNA is short. This ratio discriminates between them.

## Choosing a Workflow

| Feature | Nextflow | Snakemake |
|---|---|---|
| Language | Groovy/DSL2 | Python |
| Resumability | Built-in (`-resume`) | Built-in (`--rerun-incomplete`) |
| Container support | Native | Via `--use-singularity` |
| HPC executor | LSF / SLURM | LSF / SLURM |
| Parallelism | Automatic via channels | Checkpoint-based |

Both implementations produce identical outputs.

## Running outside Mount Sinai

The pipeline was developed on an LSF cluster with environment modules, but
neither is required.

**Tools.** Every `module load` is now conditional: if the tool is already on
`PATH` (container, conda env, system install) nothing is loaded, and if no
module system exists the call is skipped. A missing tool is reported by name
before any work starts, rather than failing mid-pipeline. Provide the tools any
way you like — see `nextflow/environment.yml`, or `env/myenv.explicit.lock.txt`
for the exact package set used in the paper.

**Scheduler.** Resource requests are declared with portable `cpus`/`memory`/`time`
directives; scheduler-specific syntax lives only in the profiles.

```bash
# Nextflow
nextflow run main.nf -profile lsf      # LSF (Minerva)
nextflow run main.nf -profile slurm    # SLURM  — set SLURM_PARTITION / SLURM_ACCOUNT
nextflow run main.nf -profile local    # single machine

# Snakemake
bash snakemake/run_snakemake.sh <run_dir>                        # LSF (default)
bash snakemake/run_snakemake.sh <run_dir> --local --extra "--cores 8"
bash snakemake/run_snakemake.sh <run_dir> --cluster-cmd 'sbatch -p {cluster.queue} -c {threads} -t {cluster.time}'
```

> **Changed default:** Nextflow used to submit to LSF even with no `-profile`.
> The executor now comes from the profile, so `-profile lsf` is required on
> Minerva; with no profile Nextflow runs locally.

The SLURM profile is **untested** — no SLURM system was available to validate it.
Treat it as a starting point and check the generated `.command.run`.

**Configuration.** Copy `snakemake/config.template.yaml`, fill in the database
paths, and pass it with `--configfile`. Required values are validated at startup
and all missing ones are reported together.

**Hardware.** KrakenUniq memory-maps the whole database: the primary MicrobialDB
step requests **180 GB RAM**. This is the main hardware barrier to running the
pipeline; the supplemental DB step needs 60 GB and everything else ≤30 GB.

## Requirements

### Software

| Tool | Version | Install |
|---|---|---|
| samtools | ≥1.18 | conda / module |
| minimap2 | ≥2.24 | conda / module |
| seqkit | ≥0.10.1 | conda / module |
| blastn (BLAST+) | ≥2.13.0 | conda / module |
| krakenuniq | 1.0.4 | conda |
| taxonkit | 0.14.1 | standalone binary (not in the conda env) |
| python | **3.7.12** | conda |
| pandas | **1.3.5** | conda |
| Nextflow | ≥22.10 | [nextflow.io](https://www.nextflow.io/) |
| Snakemake | 7.24.0 | conda |

> **Pin python and pandas.** The published results were produced with python
> 3.7.12 / pandas 1.3.5. pandas 2.x changes groupby and median semantics, both of
> which the MLA metric depends on, so a "newer is fine" environment will silently
> produce different numbers. Exact versions for all 341 packages are in
> `env/myenv.explicit.lock.txt`.

### Databases

| Database | Size | Notes |
|---|---|---|
| KrakenUniq microbial DB (primary) | 544 GB | See [KrakenUniq docs](https://github.com/fbreitwieser/krakenuniq) |
| KrakenUniq supplemental DB (viral+fungal) | 130 GB | Optional; build with `db/build_supplemental_db_v2.sh` |
| NCBI nt (BLAST) | 534 GB (179 volumes, 2024-08-31) | `update_blastdb.pl nt` |
| T2T CHM13v2 minimap2 index | 7.1 GB (`.mmi`) | 28 GB for the full reference directory |

Sizes above are measured from the Mount Sinai installation, not estimates. Total
is **~1.2 TB**, which is why these are shared in place rather than copied — see
[HANDOVER.md](HANDOVER.md).

#### Building the supplemental database

The supplemental database covers viral (all complete genomes) and fungal (complete genome + chromosome, full genome representation) assemblies from NCBI RefSeq — filling coverage gaps in the primary MicrobialDB.

```bash
# Edit paths at the top of the script, then submit to LSF:
bsub < db/build_supplemental_db_v2.sh
```

The script uses `krakenuniq-download` for automatic assembly filtering and builds the database with `krakenuniq-build`. Estimated build time: 2–4 hours on 8 threads.

## Installation

### Option 1: Conda (recommended for HPC)

```bash
conda env create -f nextflow/environment.yml
conda activate longread-microbiome
```

Then load remaining tools via environment modules:
```bash
module load samtools
module load minimap2
module load seqkit
module load blast
```

### Option 2: Singularity (most portable)

Build the container (requires internet access):
```bash
singularity build --remote nextflow/pipeline.sif nextflow/pipeline.def
```

All tools are bundled — no additional installs needed.

---

## Nextflow Usage

```bash
cd nextflow/

# LSF cluster with Singularity (recommended)
nextflow run main.nf \
    -profile lsf,singularity \
    --bam_dir /path/to/bams \
    --kraken_db /path/to/krakenuniq_db \
    --blastdb /path/to/nt \
    --t2t_ref /path/to/chm13v2.0.mmi

# LSF cluster with conda/modules
nextflow run main.nf \
    -profile lsf \
    --bam_dir /path/to/bams \
    --kraken_db /path/to/krakenuniq_db \
    --blastdb /path/to/nt \
    --t2t_ref /path/to/chm13v2.0.mmi

# Specify samples explicitly
nextflow run main.nf -profile lsf \
    --bam_dir /path/to/bams \
    --samples sample1,sample2,sample3

# Resume after failure
nextflow run main.nf -profile lsf --bam_dir /path/to/bams -resume
```

### Nextflow Parameters

| Parameter | Description | Default |
|---|---|---|
| `--bam_dir` | Directory containing input BAM files | `.` |
| `--samples` | Comma-separated sample list | auto-discover `*.bam` |
| `--resume_from_fasta` | Start from existing per-sample `*.after_t2t.unmapped.fasta.gz` and mapped-human stats files | `false` |
| `--fasta_dir` | Directory containing `SAMPLE/SAMPLE.after_t2t.unmapped.fasta.gz` for FASTA-resume mode | `.` |
| `--stats_dir` | Directory containing `SAMPLE/SAMPLE.bam.mapped_human_reads_only.stats.txt.gz` for FASTA-resume mode | `.` |
| `--kraken_db` | Path to primary KrakenUniq database | required |
| `--kraken_db_suppl` | Path to supplemental KrakenUniq database | required if `use_suppl_db=true` |
| `--use_suppl_db` | Run supplemental DB and merge results | `true` |
| `--blastdb` | Path to NCBI nt BLAST database | required |
| `--t2t_ref` | Path to T2T CHM13v2 minimap2 index | required |
| `--outdir` | Output directory | `results/` |
| `--scratch_dir` | Scratch space for BLAST temp files | `/tmp` |
| `--blast_split_nseq` | Reads per BLAST chunk | `5000` |

To skip the supplemental DB:
```bash
nextflow run main.nf -profile lsf --use_suppl_db false ...
```

To resume from existing unmapped FASTA files:
```bash
nextflow run main.nf -profile lsf \
    --resume_from_fasta true \
    --fasta_dir /path/to/sample_dirs \
    --stats_dir /path/to/sample_dirs \
    --outdir results \
    -resume
```

---

## Snakemake Usage

```bash
cd snakemake/

# Edit config.yaml to set database paths and samples
# Then run on LSF:
snakemake --snakefile snakefile \
    --cluster "bsub -P {config[lsf][project]} -q {config[lsf][queue_default]} -n {threads} -R span[hosts=1] -R rusage[mem={resources.mem_mb}]" \
    --jobs 200 \
    --configfile config.yaml

# Resume after failure
snakemake --snakefile snakefile --rerun-incomplete --configfile config.yaml
```

### Snakemake Configuration (`config.yaml`)

```yaml
samples: []                  # leave empty to auto-discover *.bam

scriptsdir: "scripts"
kraken_db: "/path/to/krakenuniq_db"
kraken_db_suppl: "/path/to/supplemental_db"  # required if use_suppl_db: true
use_suppl_db: true           # set false to skip supplemental DB
blastdb: "/path/to/nt"
t2t_ref: "/path/to/chm13v2.0.mmi"

threads:
  unmapped_analysis: 10
  krakenuniq: 4
  megablast: 16

resources:
  krakenuniq_mem_mb: 180000
  megablast_mem_mb: 30000
```

To skip the supplemental DB via CLI:
```bash
snakemake --snakefile snakefile --config use_suppl_db=false ...
```

---

## Output

Per sample:

| File | Description |
|---|---|
| `<sample>.after_t2t.unmapped.fasta.gz` | Unmapped reads after T2T filtering |
| `<sample>.krakenuniq.report` | KrakenUniq classification report |
| `<sample>.krakenuniq.info_collection.flt` | Filtered KrakenUniq hits |
| `<sample>.blast.txt` | Merged megaBLAST results |
| `<sample>.blast.microbiome.txt` | BLAST-validated microbial hits |
| `<sample>.microbiome.sum_by_length_per_genus.txt` | Genus-level abundance |
| `<sample>.median_l_adj.txt` | Fragment length-adjusted microbiome metric |

## Repository Structure

```
long_reads_microbiome_identification/
├── README.md
├── HANDOVER.md                       # running this without the original account
├── check_env.sh                      # preflight: is every DB/tool readable by you?
├── lib/
│   └── hpc_modules.sh                # conditional `module load` (no-op off Lmod)
├── env/
│   ├── myenv.explicit.lock.txt       # exact 341-package conda lock (paper env)
│   └── myenv.full_versions.txt
├── db/
│   └── build_supplemental_db_v2.sh   # build viral+fungal supplemental KrakenUniq DB
├── nextflow/
│   ├── main.nf
│   ├── nextflow.config
│   ├── environment.yml
│   ├── pipeline.def
│   ├── modules/
│   │   ├── unmapped.nf
│   │   ├── krakenuniq.nf
│   │   ├── krakenuniq_suppl.nf       # supplemental DB classification
│   │   ├── kraken_process.nf
│   │   ├── kraken_process_suppl.nf   # supplemental DB post-processing
│   │   ├── merge_kraken.nf           # merge primary + supplemental hits
│   │   ├── blast.nf
│   │   ├── blast_process.nf
│   │   ├── samtools_stats.nf
│   │   └── median_length_adj.nf
│   └── scripts/
│       ├── unmapped_analysis.sh
│       ├── get_unmapped.sh
│       ├── long_read.mm2.no_sort.sh
│       ├── krakenuniq.single.sh
│       ├── post_kraken_filter.sh
│       ├── megablast.sh
│       ├── median_length_adj.py
│       └── blast_result_process.mt.py
└── snakemake/
    ├── snakefile
    ├── config.yaml                    # Minerva paths
    ├── config.template.yaml           # start here at any other site
    ├── run_snakemake.sh
    ├── rules/
    │   ├── unmapped.smk
    │   ├── krakenuniq.smk
    │   ├── krakenuniq_suppl.smk       # supplemental DB classification
    │   ├── kraken_process.smk
    │   ├── kraken_process_suppl.smk   # supplemental DB post-processing
    │   ├── merge_kraken.smk           # merge primary + supplemental hits
    │   ├── blast.smk
    │   ├── blast_process.smk
    │   ├── samtools_stats.smk
    │   └── median_length_adj.smk
    └── scripts/
        ├── unmapped_analysis.sh
        ├── get_unmapped.sh
        ├── long_read.mm2.no_sort.sh
        ├── krakenuniq.single.sh
        ├── post_kraken_filter.sh
        ├── megablast.sh
        ├── median_length_adj.py
        └── blast_result_process.mt.py
```

## Citation

If you use this pipeline, please cite:

Zhang et al. (2026). *Critical assessment of intratumor and low-biomass microbiome using long-read sequencing.* bioRxiv. https://www.biorxiv.org/content/10.64898/2026.02.02.703393v2
