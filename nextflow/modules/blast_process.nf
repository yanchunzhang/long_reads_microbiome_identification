// modules/blast_process.nf
// Two post-BLAST processing steps:
//   PROCESS_BLAST          – multi-threaded BLAST hit processing (Python)
//   ANNOTATE_BLAST_LENGTHS – join KrakenUniq read-length info, reformat with
//                            taxonkit, and filter to microbial kingdom hits
//
// Mirrors Snakemake rules: process_blast, annotate_blast_lengths


// ── 1. Process raw BLAST output ───────────────────────────────────────────
process PROCESS_BLAST {
    tag "${sample}"

    input:
    tuple val(sample), path(blast)

    output:
    tuple val(sample), path("${sample}.blast.processed.txt")

    script:
    """
    python ${params.scriptsdir}/blast_result_process.mt.py \\
        --input   ${blast} \\
        --output  ${sample}.blast.processed.txt \\
        --threads ${task.cpus}
    """
}


// ── 2. Annotate hits with read-length info and filter to microbiome ───────
// Input channel carries: [sample, info_flt, processed_txt]
// (joined from KRAKEN_PROCESS.out.info + PROCESS_BLAST.out in main.nf)
//
// Emits two named output channels:
//   .add_length  → [sample, *.blast.processed.add_length.txt]
//   .microbiome  → [sample, *.blast.microbiome.txt]  (used downstream)

process ANNOTATE_BLAST_LENGTHS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(info), path(processed)

    output:
    tuple val(sample), path("${sample}.blast.processed.add_length.txt"), emit: add_length
    tuple val(sample), path("${sample}.blast.microbiome.txt"),           emit: microbiome

    script:
    // Note: \$ is required inside Nextflow """ blocks to pass a literal $
    // to the shell (Groovy interpolates ${...} but leaves \$ as $).
    """
    set -o pipefail

    awk 'NR==FNR {a[\$3]=\$8; next} (\$1 in a) {print \$0"\\t"a[\$1]}' \\
        ${info} ${processed} | \\
    awk '{print \$0, \$4/\$5}' | \\
    sed 's/ /\\t/g' | \\
    taxonkit lca -i 3 -s ';' -U -D | \\
    taxonkit reformat -I 7 -F -P | \\
    awk 'BEGIN {FS=OFS="\\t"} {lineage=\$8; \$7=lineage; NF=7; print}' | \\
    sed 's/ /_/g' | \\
    sort -k7,7 -k3,3 > ${sample}.blast.processed.add_length.txt

    awk '\$6>0.5 && !/k__unclass/ && !/g__unclass/ && /k__/ && \\
         !/k__Metazoa/ && \\
         (!/k__Euka/ || /(p__Ascomycota|p__Basidiomycota|p__Mucoromycota|p__Chytridiomycota)/)' \\
        ${sample}.blast.processed.add_length.txt > ${sample}.blast.microbiome.txt
    """
}
