#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Module imports ────────────────────────────────────────────────────────
include { UNMAPPED_ANALYSIS }       from './modules/unmapped'
include { MAPPED_HUMAN_STATS }      from './modules/samtools_stats'
include { KRAKENUNIQ }              from './modules/krakenuniq'
include { KRAKENUNIQ_SUPPL }        from './modules/krakenuniq_suppl'
include { KRAKEN_PROCESS }          from './modules/kraken_process'
include { KRAKEN_PROCESS_SUPPL }    from './modules/kraken_process_suppl'
include { MERGE_KRAKEN }            from './modules/merge_kraken'
include { SPLIT_BLAST_QUERY;
          MEGABLAST_CHUNK;
          MERGE_BLAST_CHUNKS }      from './modules/blast'
include { PROCESS_BLAST;
          ANNOTATE_BLAST_LENGTHS }  from './modules/blast_process'
include { MEDIAN_LENGTH_ADJ }       from './modules/median_length_adj'

// ── Workflow ──────────────────────────────────────────────────────────────
workflow {

    // Build sample list when explicitly provided.
    // Handles both a config list (["A1","A2"]) and a CLI comma-separated
    // string (--samples A1,A2,A3) by tokenising on commas when needed.
    def sample_list = []
    if (params.samples) {
        sample_list = (params.samples instanceof List)
            ? params.samples.collect { it.toString().trim() }
            : params.samples.toString().tokenize(',').collect { it.trim() }
    }

    if (params.resume_from_fasta.toString().toLowerCase() == 'true') {
        // FASTA-resume mode mirrors Snakemake's ability to start from existing
        // per-sample unmapped FASTA files. Precomputed human-mapped stats are
        // still required for the final median-length adjustment step.
        if (sample_list) {
            ch_unmapped_fasta = Channel
                .from(sample_list)
                .map { sample -> tuple(sample, file("${params.fasta_dir}/${sample}/${sample}.after_t2t.unmapped.fasta.gz")) }

            ch_stats = Channel
                .from(sample_list)
                .map { sample -> tuple(sample, file("${params.stats_dir}/${sample}/${sample}.bam.mapped_human_reads_only.stats.txt.gz")) }
        } else {
            ch_unmapped_fasta = Channel
                .fromPath("${params.fasta_dir}/*/*.after_t2t.unmapped.fasta.gz")
                .map { fa -> tuple(fa.parent.name, fa) }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .flatMap()

            ch_stats = Channel
                .fromPath("${params.stats_dir}/*/*.bam.mapped_human_reads_only.stats.txt.gz")
                .map { stats -> tuple(stats.parent.name, stats) }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .flatMap()
        }
    } else {
        // Build sample → BAM channel
        // If params.samples is non-empty, use it; otherwise auto-discover */*.bam.
        if (sample_list) {
            ch_bam = Channel
                .from(sample_list)
                .map { sample -> tuple(sample, file("${params.bam_dir}/${sample}/${sample}.bam")) }
        } else {
            ch_bam = Channel
                .fromPath("${params.bam_dir}/*/*.bam")
                .filter { bam -> !bam.name.contains("unmapp") && !bam.name.contains("mm2") }
                .map    { bam -> tuple(bam.parent.name, bam) }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .flatMap()
        }

        // Step 1 (parallel): extract unmapped reads + compute human-mapping stats
        UNMAPPED_ANALYSIS(ch_bam)
        MAPPED_HUMAN_STATS(ch_bam)

        ch_unmapped_fasta = UNMAPPED_ANALYSIS.out
        ch_stats = MAPPED_HUMAN_STATS.out
    }

    // Step 2: KrakenUniq — primary DB always runs; supplemental is optional
    KRAKENUNIQ(ch_unmapped_fasta)
    KRAKEN_PROCESS(KRAKENUNIQ.out)

    if (params.use_suppl_db.toString().toLowerCase() != 'false') {
        // Step 2b: Supplemental DB in parallel
        KRAKENUNIQ_SUPPL(ch_unmapped_fasta)
        KRAKEN_PROCESS_SUPPL(KRAKENUNIQ_SUPPL.out)

        // Step 3: Merge primary + supplemental (deduplicate by read ID)
        ch_merge = KRAKEN_PROCESS.out.fasta
            .join(KRAKEN_PROCESS_SUPPL.out.fasta)
            .join(KRAKEN_PROCESS.out.info)
            .join(KRAKEN_PROCESS_SUPPL.out.info)
            .join(KRAKENUNIQ.out.map { sample, krakenuniq, classified, report -> tuple(sample, krakenuniq) })
        MERGE_KRAKEN(ch_merge)

        ch_blast_fasta   = MERGE_KRAKEN.out.fasta
        ch_annotate_info = MERGE_KRAKEN.out.info
    } else {
        ch_blast_fasta   = KRAKEN_PROCESS.out.fasta
        ch_annotate_info = KRAKEN_PROCESS.out.info
    }

    // Step 4: Split FASTA → run megaBLAST per chunk (in parallel) → merge
    SPLIT_BLAST_QUERY(ch_blast_fasta)

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
    ch_annotate = ch_annotate_info
        .join(PROCESS_BLAST.out)

    ANNOTATE_BLAST_LENGTHS(ch_annotate)

    // Step 7: Median-length-adjusted microbiome abundance
    ch_final = ch_stats
        .join(ANNOTATE_BLAST_LENGTHS.out.microbiome)

    MEDIAN_LENGTH_ADJ(ch_final)
}
