# Supplemental KrakenUniq Database

## Background

The primary database (`kuniq_microbialdb_minus_kdb.20230808`, 535 GB) covers bacteria
and a subset of eukaryotic pathogens from EuPathDB. A detailed audit of its
`database.kdb.counts` and `seqid2taxid.map` revealed two significant coverage gaps:

### What the primary database contains

Confirmed by direct audit of `database.kdb.counts` and `seqid2taxid.map`
(2025-05-14). The primary DB is the Microbial2023 database (Salzberg et al.):

| Kingdom | Status | Notes |
|---|---|---|
| Bacteria | ✓ Good | Core of the DB; 34,452 genomes (16,150 species) |
| Archaea | ✓ Sparse | 422 genomes; a few environmental genera |
| Viruses | ⚠ Gaps | ~12,910 species total, but spot-check confirmed missing: Influenza A/B, HCV, HBV, HSV-2, JC polyomavirus, Human adenovirus |
| Fungi — clinical | ✓ Present | Candida albicans, C. auris, C. glabrata, Cryptococcus neoformans, Pneumocystis jirovecii — **but were being discarded by the Eukaryota filter (now fixed)** |
| Fungi — general | ⚠ EuPathDB only | Dominated by plant pathogens (Hemileia, Puccinia) and microsporidia; limited broader fungal coverage |
| EuPathDB protozoa | ✓ Present | Toxoplasma, Trypanosoma, Plasmodium, Acanthamoeba, etc. |

### What the supplemental database adds

- **RefSeq viral genomes (Complete + Chromosome level)** — fills confirmed gaps: Influenza A/B, HCV, HBV, HSV-2, JC polyomavirus, and updates post-2023 strains
- **RefSeq fungal genomes (Complete + Chromosome level)** — broader coverage beyond EuPathDB scope

---

## Genome selection criteria

Consistent with Salzberg et al. (*Sci Transl Med* 2025) and Dohlman et al.
PathSeq-T2T (*Cell* 2026), both of which use complete/chromosome-level
RefSeq assemblies as their reference standard.

|  | Viruses | Fungi |
|---|---|---|
| Source | NCBI RefSeq | NCBI RefSeq |
| `assembly_level` | **Complete Genome only** | Complete Genome or Chromosome |
| `genome_rep` | — (not filtered) | Full |
| Min contig length | None (complete genomes are full-length by definition) | ≥ 5,000 bp |

**Viral filter rationale:** Salzberg uses "all RefSeq viral *complete genomes*" —
the strictest standard. Complete viral genomes are inherently single contiguous
sequences so chromosome-level is not applicable.

**Salzberg vs PathSeq-T2T on fungi:** Salzberg used *all* 572 RefSeq fungal
assemblies (no assembly-level filter) for their Fungi_RefSeq database, while
PathSeq-T2T (PlusPF) uses complete/chromosome-level only. We follow
PathSeq-T2T here to reduce false positives from fragmented scaffolds.

### What is NOT included (and why)
- **Scaffold/Contig-level assemblies**: higher contamination rate; Contig N50
  typically <50 kbp for fungi.
- **GenBank-only assemblies**: not subject to NCBI FCS-GX contamination screening.
- **MAGs**: ~1–5% cross-species contamination even in high-quality bins.
- **Bacteria / Archaea**: already comprehensively covered in the primary DB.

---

## Build history

### v1 — `kuniq_supplemental_vf_db` (2026-05-15, completed)

**Script:** `db/build_supplemental_db.sh`  
**Approach:** Manual download loop using `wget` against NCBI assembly summaries, then
`krakenuniq-build --add-to-library` per FASTA.

**Issues encountered:**
- `krakenuniq-build --add-to-library` exits with code 255 even on success (KrakenUniq
  v1.0.4 Perl bug) — had to use `|| true` to prevent job abort.
- `krakenuniq-build --build` (step 4) failed with `JELLYFISH_BIN: unbound variable`
  because its internal `set -u` was active — fixed with
  `export JELLYFISH_BIN="${JELLYFISH_BIN:-jellyfish}"`.
- `count_unique` exceeded CPU count because LSF spread 8 threads across multiple nodes
  — fixed by adding `#BSUB -R "span[hosts=1]"`.
- `seqid2taxid.map` was empty after the crash: the `--add-to-library` step never wrote
  `.map` files before exiting 255. Manually reconstructed from `|kraken:taxid|` FASTA
  headers using:
  ```bash
  awk 'match($0, /\|kraken:taxid\|([0-9]+)/, arr) { id=substr($0,2); gsub(/ .*/,"",id); print id"\t"arr[1] }' library/added/*.fna > seqid2taxid.map
  ```
  (First attempt used `awk -F'|kraken:taxid|'` which failed because `|` is regex OR.)

**Result:** 39 GB database.kdb, 205K database.kdb.counts, 14,302 unique taxids.
`seqid2taxid.map` is 783K (manually reconstructed — covers only viral FASTA headers,
not all contigs).

---

### v2 — `kuniq_supplemental_vf_db_v2` (2026-05-15, completed) ← **current**

**Script:** `db/build_supplemental_db_v2.sh`  
**Approach:** Uses `krakenuniq-download` natively, which handles assembly filtering,
downloading, and `seqid2taxid.map` generation in one step. Then `krakenuniq-build --build`.

```bash
krakenuniq-download --db "$SUPDB" --threads 8 refseq/viral/Complete_Genome
krakenuniq-download --db "$SUPDB" --threads 8 --min-seq-len 5000 'refseq/fungi/Complete_Genome/genome_rep=Full'
krakenuniq-download --db "$SUPDB" --threads 8 --min-seq-len 5000 'refseq/fungi/Chromosome/genome_rep=Full'
export JELLYFISH_BIN="${JELLYFISH_BIN:-jellyfish}"
krakenuniq-build --build --db "$SUPDB" --threads 8 --kmer-len 31 --minimizer-len 15
```

No crashes. `seqid2taxid.map` built natively and contains all contig-level mappings.

**Result:** 39 GB database.kdb, 205K database.kdb.counts, **14,305 unique taxids**.
`seqid2taxid.map` is 1.4M (native, complete).

**Why v2 is preferred over v1:** identical k-mer content, but the native `seqid2taxid.map`
is complete (all contigs mapped) vs. the manually reconstructed v1 map which only
captured viral header-level taxids. This ensures correct LCA assignment for all sequences.

---

## Current database path
```
/sc/arion/projects/schzrnas/zhangy40/softwares/kuniq_supplemental_vf_db_v2
```

## Build script (recommended)

```bash
bsub < db/build_supplemental_db_v2.sh
```

Estimated resources: 64 GB RAM, 8 cores, ~4 h wall time.

---

## Pipeline integration

The supplemental database is run **in addition to** the primary database, not
as a replacement. Reads classified by either database are retained.

### Nextflow — add to `nextflow.config`

```groovy
params {
    kraken_db       = "/sc/arion/projects/schzrnas/zhangy40/softwares/kuniq_microbialdb_minus_kdb.20230808"
    kraken_db_suppl = "/sc/arion/projects/schzrnas/zhangy40/softwares/kuniq_supplemental_vf_db_v2"
    use_suppl_db    = true   // set false to skip
}
```

### Post-KrakenUniq filter — `post_kraken_filter.sh` ✓ updated

Both `nextflow/scripts/post_kraken_filter.sh` and
`snakemake/scripts/post_kraken_filter.sh` have been updated (line 9).

```bash
# Before — dropped ALL Eukaryota, silently discarding EuPathDB pathogens,
# clinical fungi (Candida, Cryptococcus, Pneumocystis), and viruses that
# resolve to a eukaryotic host in the taxonomy
awk '!/Eukaryota/ && /d__/ && !/synthetic/'

# After — retains eukaryotic pathogens; excludes only human reads and
# synthetic/spike-in sequences
awk '/d__/ && !/synthetic/ && !/Homo sapiens/'
```

Impact: this change alone unlocks detection of the following organisms that
were already indexed in the primary database but were being discarded:
Candida albicans, C. auris, C. glabrata, Cryptococcus neoformans,
Pneumocystis jirovecii, Toxoplasma gondii, Trypanosoma cruzi, CMV, EBV, HSV-1,
SARS-CoV-2, and all EuPathDB protozoa.

---

## Insights from literature review (2025-05-14)

### Salzberg et al., *Sci Transl Med* 2025
- Confirmed our primary DB **is** Microbial2023 (same 535 GB, same contents)
- Two-pass host depletion: GRCh38 → CHM13 T2T — our pipeline already does this ✓
- Separate fungi database strategy: second KrakenUniq pass on reads that were
  *unclassified OR classified as fungi* (taxid 4751); used all 572 RefSeq fungi
- Used KrakenUniq default parameters; no per-read unique k-mer threshold

### Dohlman et al. (PathSeq-T2T), *Cell* 2026
- Multi-step host filter adds: MHC/HLA loci, breakpoint junctions, UniVec,
  Gencode transcripts before T2T alignment — primarily beneficial for short reads
- Three classifiers: Kraken2 (`--confidence 0.15`) + MetaPhlAn4 + Sylph
- Database: PlusPF (bacteria/archaea/viruses/fungi/protozoa, complete +
  chromosome level) — this is the basis for our assembly-level filter choice
- RPM normalization: divide by **total primary reads** (human + microbial),
  not just unmapped reads
- Contamination background: ~1 RPM from cross-contamination in sterile samples;
  common lab contaminants to watch: *Cutibacterium acnes*, *Staphylococcus
  epidermidis*, *Malassezia restricta*, *Pseudomonas*, *Sphingomonas*,
  *Bradyrhizobium*

### Changes applied based on these papers
| Change | File | Reason |
|---|---|---|
| Remove `!/Eukaryota/` filter | `post_kraken_filter.sh` (both) | Clinical fungi and EuPathDB protozoa already in primary DB were being discarded |
| Viral DB: **Complete Genome only** | `build_supplemental_db.sh` | Matches Salzberg exactly; previously used all assembly levels |

---

## Key references

- Walker MA et al. (2018) PathSeq: a customizable computational tool for the
  identification of known and novel pathogens. *Cell Systems* 7(5).
  https://doi.org/10.1016/j.cels.2018.09.016
- Breitwieser FP et al. (2018) KrakenUniq: confident and fast metagenomics
  classification using unique k-mer counts. *Genome Biology* 19, 198.
  https://doi.org/10.1186/s13059-018-1568-0
- NCBI RefSeq selection criteria:
  https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/
- NCBI FCS-GX contamination screening:
  https://github.com/ncbi/fcs
