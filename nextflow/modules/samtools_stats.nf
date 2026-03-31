// modules/samtools_stats.nf
// Computes samtools stats for human-mapped reads only.
// Mirrors Snakemake rule: mapped_human_stats

process MAPPED_HUMAN_STATS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.bam.mapped_human_reads_only.stats.txt")

    script:
    """
    set -euo pipefail

    sh ${params.scriptsdir}/samtools_stats.mapped_human_reads_only.sh \\
        ${bam} \\
        ${task.cpus}
    """
}
