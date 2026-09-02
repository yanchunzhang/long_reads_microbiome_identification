# Running the test dataset

A smoke test that exercises all 13 pipeline rules on a small BAM. Run this once
on a new account, before committing a real sample to the queue.

Contents of this directory:

| File | |
|---|---|
| `CO4_T1.test.bam` | 19,899 reads — 19,836 mapped human + **63 unmapped** |
| `CO4_T1.test.bam.bai` | index |
| `create_test_bam.sh` | how the BAM was built (see the warning at the end) |

Only the **63 unmapped** reads feed the microbiome branch; the mapped reads
exercise the human-stats rule and provide the host median length that MLA is
calculated against. Read lengths of the unmapped set span 45–12,896 bp
(median 3,453), so the long-read behaviour is genuinely tested.

---

## Before you start

```bash
bash check_env.sh
```

Expect `passed 23  warnings 0  failures 0`. If databases fail here, stop and fix
that first — see [HANDOVER.md](../HANDOVER.md). Everything below assumes the
preflight is green.

## Run it

Pick the repo location you have access to:

```bash
# Fang lab (fangg03a) — self-contained mirror, no schzrnas access needed
P=/sc/arion/projects/fangg03a/zhangy40/long_read_microbiome/pipeline

# or, with schzrnas access — the canonical repo
P=/sc/arion/projects/schzrnas/zhangy40/scripts/long_reads_microbiome_identification/long_reads_microbiome_identification
```

Then:

```bash
R=/sc/arion/scratch/$USER/pipeline_test_run

mkdir -p $R/CO4_T1
cp $P/test/CO4_T1.test.bam     $R/CO4_T1/CO4_T1.bam
cp $P/test/CO4_T1.test.bam.bai $R/CO4_T1/CO4_T1.bam.bai

bash $P/snakemake/run_snakemake.sh $R -n  -A acc_fangg03a   # dry run: expect 13 rules
bash $P/snakemake/run_snakemake.sh $R -j 20 -A acc_fangg03a # real run
```

## Four things that catch people out

**1. The file must be renamed on the way in.** Samples are discovered as
`<run_dir>/<sample>/<sample>.bam`. The test ships as `CO4_T1.test.bam` and must
land as `CO4_T1/CO4_T1.bam`. Copy it under its own name and the sample becomes
`CO4_T1.test`, the expected `.bam` never matches, and Snakemake reports *nothing
to do* rather than an error — the most confusing possible failure.

**2. `-A` must match your LSF billing account.** The default is `acc_schzrnas`.
Fang lab members need `-A acc_fangg03a` (or `export LSF_ACCOUNT=acc_fangg03a`).
The Unix group and the LSF account are different things: `check_env.sh` will
pass and then every `bsub` is rejected at submission.

**3. Run inside `tmux` or `screen`.** The launcher stays alive managing the DAG
for the whole run; if your SSH session drops, the run dies. To recover:
`bash $P/snakemake/run_snakemake.sh $R -u` unlocks the directory, then re-run.

**4. Write to scratch, not into the pipeline directory.** Outputs go back into
the run directory. Use `/sc/arion/scratch/$USER/...`.

## How long it takes

Roughly an hour, and **it is not proportional to the test size**. Runtime is
dominated by KrakenUniq loading the 544 GB primary database into a 180 GB
high-memory node. A 63-read input costs nearly what a real sample costs, so
don't read a slow run as a problem.

Rules run mostly one at a time for a single sample; `bjobs` will usually show
one job. The megablast step fans out into chunks.

## What you should see

Progress in order: `mapped_human_stats` → `unmapped_analysis` → `krakenuniq` +
`krakenuniq_suppl` → `kraken_process*` → `merge_kraken` →
`split_blast_query_fasta` → `megablast_chunk` → `merge_blast_chunks` →
`process_blast` → `annotate_blast_lengths` → `median_length_adj`.

Key outputs in `$R/CO4_T1/`:

```
CO4_T1.after_t2t.unmapped.fasta.gz          reads surviving both host filters
CO4_T1.krakenuniq / .suppl.krakenuniq       primary and supplemental classification
CO4_T1.merged.krakenuniq.microbiome.fasta   merged, human/synthetic removed
CO4_T1.blast.microbiome.txt                 validated microbial reads
CO4_T1.microbiome.sum_by_length_per_genus.txt
CO4_T1.median_l_adj.txt                     the MLA values
CO4_T1.human_median_length.tsv              host median (MLA denominator)
```

Compare against `test/expected/` (see below) to confirm not just that it ran,
but that it produced the right answer.

---

## WARNING: this BAM cannot be regenerated

`create_test_bam.sh` reads from
`intratumor_bacteria/snakemake/CO4_T1.bam` and `CO4_T1.blast.microbiome.txt`.
**That directory was deleted on 2026-08-07** and both source files are gone. The
script is kept for provenance — it documents how the BAM was built — but it will
not run.

`CO4_T1.test.bam` is therefore irreplaceable. Copies:

- this repo (and the fangg03a mirror)
- `onedrive:BACTERIA_IN_TUMOR/pipeline_test_data/` (verified by `rclone check`)

It is `.gitignore`d (138 MB exceeds GitHub's 100 MB file limit), so it does not
ship with the public repository. Do not delete it from any of the above without
confirming another copy exists.
