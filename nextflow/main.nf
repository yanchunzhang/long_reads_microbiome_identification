#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Module imports ────────────────────────────────────────────────────────
include { UNMAPPED_ANALYSIS }       from './modules/unmapped'
include { MAPPED_HUMAN_STATS }      from './modules/samtools_stats'
include { KRAKENUNIQ }              from './modules/krakenuniq'
include { KRAKEN_PROCESS }          from './modules/kraken_process'
include { SPLIT_BLAST_QUERY;
          MEGABLAST_CHUNK;
          MERGE_BLAST_CHUNKS }      from './modules/blast'
include { PROCESS_BLAST;
          ANNOTATE_BLAST_LENGTHS }  from './modules/blast_process'
include { MEDIAN_LENGTH_ADJ }       from './modules/median_length_adj'

// ── Workflow ──────────────────────────────────────────────────────────────
workflow {

    // Build sample → BAM channel
    // If params.samples is non-empty, use it; otherwise auto-discover *.bam
    // Handles both a config list (["A1","A2"]) and a CLI comma-separated
    // string (--samples A1,A2,A3) by tokenising on commas when needed.
    if (params.samples) {
        def sample_list = (params.samples instanceof List)
            ? params.samples
            : params.samples.toString().tokenize(',')

        ch_bam = Channel
            .from(sample_list)
            .map { sample -> tuple(sample.trim(), file("${params.bam_dir}/${sample.trim()}.bam")) }
    } else {
        ch_bam = Channel
            .fromPath("${params.bam_dir}/*.bam")
            .filter { bam -> !bam.name.contains("unmapp") && !bam.name.contains("mm2") }
            .map    { bam -> tuple(bam.simpleName, bam) }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .flatMap()
    }

    // Step 1 (parallel): extract unmapped reads + compute human-mapping stats
    UNMAPPED_ANALYSIS(ch_bam)
    MAPPED_HUMAN_STATS(ch_bam)

    // Step 2: KrakenUniq taxonomic classification on unmapped FASTA
    KRAKENUNIQ(UNMAPPED_ANALYSIS.out)

    // Step 3: Filter / post-process KrakenUniq output
    KRAKEN_PROCESS(KRAKENUNIQ.out)

    // Step 4: Split microbiome FASTA → run megaBLAST per chunk (in parallel) → merge
    SPLIT_BLAST_QUERY(KRAKEN_PROCESS.out.fasta)

    // Record per-sample chunk counts from the split directory
    ch_counts = SPLIT_BLAST_QUERY.out
        .map { sample, split_dir ->
            def n = split_dir.listFiles().findAll { it.name.endsWith('.fa') }.size()
            tuple(sample, n)
        }

    // Expand split directory into one [sample, fa_file] tuple per chunk
    ch_chunks = SPLIT_BLAST_QUERY.out
        .flatMap { sample, split_dir ->
            split_dir.listFiles()
                .findAll { it.name.endsWith('.fa') }
                .sort()
                .collect { fa -> [sample, fa] }
        }

    MEGABLAST_CHUNK(ch_chunks)

    // Collect all chunks for a sample, then merge into a single BLAST result.
    // Use combine(by:0) instead of join() — join() consumes the right-channel
    // entry on first match so only 1-chunk samples get a groupKey. combine(by:0)
    // broadcasts ch_counts to every matching MEGABLAST_CHUNK entry per sample.
    ch_grouped = MEGABLAST_CHUNK.out
        .combine(ch_counts, by: 0)
        .map { sample, blast, n -> tuple(groupKey(sample, n), blast) }
        .groupTuple()

    MERGE_BLAST_CHUNKS(ch_grouped)

    // Step 5: Process raw BLAST hits
    PROCESS_BLAST(MERGE_BLAST_CHUNKS.out)

    // Step 6: Annotate BLAST hits with KrakenUniq length/taxonomy info
    ch_annotate = KRAKEN_PROCESS.out.info
        .join(PROCESS_BLAST.out)

    ANNOTATE_BLAST_LENGTHS(ch_annotate)

    // Step 7: Median-length-adjusted microbiome abundance
    ch_final = MAPPED_HUMAN_STATS.out
        .join(ANNOTATE_BLAST_LENGTHS.out.microbiome)

    MEDIAN_LENGTH_ADJ(ch_final)
}
