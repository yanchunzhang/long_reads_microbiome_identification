// modules/samtools_stats.nf
// Computes samtools stats for human-mapped reads only.
// Mirrors Snakemake rule: mapped_human_stats

process MAPPED_HUMAN_STATS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.bam.mapped_human_reads_only.stats.txt.gz")

    script:
    """
    set -euo pipefail

    samtools view -@ ${task.cpus} ${bam} -h -F4 | \\
        samtools stats | \\
        gzip > ${sample}.bam.mapped_human_reads_only.stats.txt.gz
    """
}
