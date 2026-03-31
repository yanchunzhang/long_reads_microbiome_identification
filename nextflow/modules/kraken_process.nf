// modules/kraken_process.nf
// Filters KrakenUniq output to microbial reads and collects read-level info.
// Mirrors Snakemake rule: kraken_process
//
// Emits two named output channels:
//   .fasta  → [sample, microbiome.fasta]   fed into BLAST splitting
//   .info   → [sample, info_collection.flt] fed into ANNOTATE_BLAST_LENGTHS

process KRAKEN_PROCESS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(krakenuniq), path(classified)

    output:
    tuple val(sample), path("${sample}.krakenuniq.microbiome.fasta"), emit: fasta
    tuple val(sample), path("${sample}.krakenuniq.info_collection.flt"),  emit: info

    script:
    // The upstream shell script operates on files named by sample prefix
    // and expects the KrakenUniq report to be present in the work directory.
    """
    sh ${params.scriptsdir}/post_kraken_filter.sh ${sample} ${params.kraken_db}
    """
}
