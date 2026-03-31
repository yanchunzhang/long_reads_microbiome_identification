// modules/unmapped.nf
// Extracts unmapped (non-human) reads from a BAM and produces a gzipped FASTA.
// Mirrors Snakemake rule: unmapped_analysis

process UNMAPPED_ANALYSIS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.after_t2t.unmapped.fasta.gz")

    script:
    """
    sh ${params.scriptsdir}/unmapped_analysis.sh \\
        ${sample} ${bam} ${task.cpus} ${params.scriptsdir} ${params.t2t_ref}
    """
}
