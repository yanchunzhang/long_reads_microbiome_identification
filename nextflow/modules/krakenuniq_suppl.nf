// modules/krakenuniq_suppl.nf
// Runs KrakenUniq on the supplemental database (viral + fungal extensions).
// Mirrors KRAKENUNIQ but uses params.kraken_db_suppl and prefixes outputs
// with ".suppl" so they coexist with the primary-DB outputs per sample.

process KRAKENUNIQ_SUPPL {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}.suppl.krakenuniq"), path("${sample}.suppl.classified_by_krakenuniq")

    script:
    """
    sh ${params.scriptsdir}/krakenuniq.single.sh \\
        ${params.kraken_db_suppl} \\
        ${fasta} \\
        ${sample}.suppl \\
        ${task.cpus}
    """
}
