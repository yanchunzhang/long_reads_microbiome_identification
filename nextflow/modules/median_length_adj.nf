// modules/median_length_adj.nf
// Computes median-read-length-adjusted microbiome abundances.
// Mirrors Snakemake rule: median_length_adj
//
// Input channel carries: [sample, stats_txt, microbiome_txt]
// (joined from MAPPED_HUMAN_STATS.out + ANNOTATE_BLAST_LENGTHS.out.microbiome)

process MEDIAN_LENGTH_ADJ {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(stats), path(microbe)

    output:
    tuple val(sample),
          path("${sample}.median_l_adj.txt"),
          path("${sample}.microbiome.sum_by_length_per_genus.txt")

    script:
    """
    python ${params.scriptsdir}/median_length_adj.py \\
        --stats   ${stats} \\
        --microbe ${microbe} \\
        --sum_out ${sample}.microbiome.sum_by_length_per_genus.txt \\
        --gt5_out ${sample}.median_l_adj.txt \\
        --sample  ${sample}
    """
}
