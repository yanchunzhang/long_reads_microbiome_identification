// modules/kraken_process_suppl.nf
// Filters supplemental-DB KrakenUniq output to microbial reads.
// Mirrors KRAKEN_PROCESS but operates on the .suppl-prefixed files.

process KRAKEN_PROCESS_SUPPL {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(krakenuniq), path(classified)

    output:
    tuple val(sample), path("${sample}.suppl.krakenuniq.microbiome.fasta"), emit: fasta
    tuple val(sample), path("${sample}.suppl.krakenuniq.info_collection.flt"),  emit: info

    script:
    """
    sh ${params.scriptsdir}/post_kraken_filter.sh ${sample}.suppl ${params.kraken_db_suppl}
    """
}
