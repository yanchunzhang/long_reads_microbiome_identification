// modules/krakenuniq.nf
// Runs KrakenUniq taxonomic classification on unmapped reads.
// Mirrors Snakemake rule: krakenuniq

process KRAKENUNIQ {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}.krakenuniq"), path("${sample}.classified_by_krakenuniq")

    script:
    """
    sh ${params.scriptsdir}/krakenuniq.single.sh \\
        ${params.kraken_db} \\
        ${fasta} \\
        ${sample} \\
        ${task.cpus}
    """
}
