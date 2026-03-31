// modules/blast.nf
// Three processes that together replace the Snakemake checkpoint pattern:
//   1. SPLIT_BLAST_QUERY  – split microbiome FASTA into fixed-size chunks
//   2. MEGABLAST_CHUNK    – run megaBLAST on a single chunk (parallelised)
//   3. MERGE_BLAST_CHUNKS – concatenate per-chunk results into one file
//
// The Snakemake checkpoint is replaced by:
//   SPLIT_BLAST_QUERY emits a directory → .flatMap in main.nf lists the
//   individual chunk files and fans them out into parallel MEGABLAST_CHUNK jobs.


// ── 1. Split FASTA into chunks ─────────────────────────────────────────────
// Mirrors Snakemake checkpoint: split_blast_query_fasta

process SPLIT_BLAST_QUERY {
    tag "${sample}"

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("split_fasta_${sample}"), emit: split_dir

    script:
    """
    set -euo pipefail

    mkdir -p split_fasta_${sample}

    seqkit split2 \\
        -s ${params.blast_split_nseq} \\
        -O split_fasta_${sample} \\
        ${fasta}

    # Rename chunks to a consistent zero-padded pattern: {sample}.part_NNN.fa
    i=0
    found=0
    for f in split_fasta_${sample}/*.fasta split_fasta_${sample}/*.fa; do
        [ -e "\$f" ] || continue
        found=1
        new=\$(printf "split_fasta_${sample}/${sample}.part_%03d.fa" "\$i")
        mv "\$f" "\$new"
        i=\$((i+1))
    done

    if [ "\$found" -eq 0 ]; then
        echo "ERROR: no split FASTA files were produced for ${sample}" >&2
        exit 1
    fi
    """
}


// ── 2. Run megaBLAST on a single chunk ────────────────────────────────────
// Mirrors Snakemake rule: megablast_chunk
// Output is marked temporary in the original; publishDir is intentionally
// omitted here – only the merged result is published.

process MEGABLAST_CHUNK {
    tag "${sample} / ${fa.simpleName}"

    input:
    tuple val(sample), path(fa)

    output:
    tuple val(sample), path("${fa.baseName}.blast.txt")

    script:
    """
    set -euo pipefail

    sh ${params.scriptsdir}/megablast.sh \\
        ${fa} \\
        ${params.blastdb} \\
        ${fa.baseName}.blast.txt \\
        ${task.cpus} blast
    """
}


// ── 3. Merge all per-chunk BLAST results ─────────────────────────────────
// Mirrors Snakemake rule: merge_blast_chunks
// Chunks are sorted by filename before concatenation to match the original
// Python-based merge logic (sorted(input.chunks)).

process MERGE_BLAST_CHUNKS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(chunks)

    output:
    tuple val(sample), path("${sample}.blast.txt")

    script:
    // Sort chunks by name in Groovy so the shell cat is deterministic
    def sorted = (chunks instanceof List ? chunks : [chunks])
        .sort { a, b -> a.name <=> b.name }
    """
    cat ${sorted.join(' ')} > ${sample}.blast.txt
    """
}
