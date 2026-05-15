# long_reads_microbiome_identification

A pipeline for detecting genuine microbiome signals in low-biomass tissues using long-read sequencing (ONT/PacBio). Available in both **Nextflow** and **Snakemake** implementations.

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

## Requirements

### Software

| Tool | Version | Install |
|---|---|---|
| samtools | ≥1.18 | conda / module |
| minimap2 | ≥2.24 | conda / module |
| seqkit | ≥0.10.1 | conda / module |
| blastn (BLAST+) | ≥2.13.0 | conda / module |
| krakenuniq | ≥1.0.4 | conda |
| taxonkit | ≥0.14.1 | conda |
| python | ≥3.8 | conda |
| pandas | ≥2.0 | conda |
| Nextflow | ≥22.10 | [nextflow.io](https://www.nextflow.io/) |
| Snakemake | ≥7.0 | conda |

### Databases

| Database | Size | Notes |
|---|---|---|
| KrakenUniq microbial DB (primary) | ~180 GB | See [KrakenUniq docs](https://github.com/fbreitwieser/krakenuniq) |
| KrakenUniq supplemental DB (viral+fungal) | ~47 GB | Optional; build with `db/build_supplemental_db_v2.sh` |
| NCBI nt (BLAST) | ~500 GB | `update_blastdb.pl nt` |
| T2T CHM13v2 reference | ~3 GB | Pre-built minimap2 index (`.mmi`) recommended |

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
    ├── config.yaml
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
