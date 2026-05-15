// modules/merge_kraken.nf
// Merges primary-DB and supplemental-DB KrakenUniq outputs before BLAST.
// Deduplicates by read ID — primary DB takes priority on any read classified
// by both databases.

process MERGE_KRAKEN {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample),
          path(primary_fasta), path(suppl_fasta),
          path(primary_info),  path(suppl_info)

    output:
    tuple val(sample), path("${sample}.merged.krakenuniq.microbiome.fasta"), emit: fasta
    tuple val(sample), path("${sample}.merged.krakenuniq.info_collection.flt"),  emit: info

    script:
    """
    # Merge FASTA — skip any read ID already seen (primary DB has priority)
    awk '/^>/ { if (seen[\$0]) { skip=1 } else { seen[\$0]=1; skip=0; print } } \\
         !/^>/ && !skip { print }' \\
        ${primary_fasta} ${suppl_fasta} > ${sample}.merged.krakenuniq.microbiome.fasta

    # Merge info collection — first occurrence wins (primary DB first)
    awk '!seen[\$1]++' ${primary_info} ${suppl_info} \\
        > ${sample}.merged.krakenuniq.info_collection.flt
    """
}
