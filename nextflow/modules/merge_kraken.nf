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
          path(primary_info),  path(suppl_info),
          path(primary_kraken)

    output:
    tuple val(sample), path("${sample}.merged.krakenuniq.microbiome.fasta"), emit: fasta
    tuple val(sample), path("${sample}.merged.krakenuniq.info_collection.flt"),  emit: info

    script:
    """
    set -euo pipefail

    # Merge FASTA — skip any read ID already seen (primary DB has priority)
    awk '/^>/ { if (seen[\$0]) { skip=1 } else { seen[\$0]=1; skip=0; print } } \\
         !/^>/ && !skip { print }' \\
        ${primary_fasta} ${suppl_fasta} > ${sample}.merged.krakenuniq.microbiome.fasta.tmp

    # Merge info collection — first occurrence wins (primary DB first)
    awk '!seen[\$1]++' ${primary_info} ${suppl_info} \\
        > ${sample}.merged.krakenuniq.info_collection.flt.tmp

    # Remove reads classified as human or synthetic by the primary DB.
    awk '\$3==9606 || \$3==32630 {print \$2}' ${primary_kraken} \\
        > ${sample}.merged.krakenuniq.microbiome.fasta.human_synthetic.txt

    if [[ -s "${sample}.merged.krakenuniq.microbiome.fasta.human_synthetic.txt" ]]; then
        seqkit grep -v -f ${sample}.merged.krakenuniq.microbiome.fasta.human_synthetic.txt \\
            ${sample}.merged.krakenuniq.microbiome.fasta.tmp \\
            > ${sample}.merged.krakenuniq.microbiome.fasta
        awk 'NR==FNR {seen[">"\$1]=1; next} !(\$1 in seen)' \\
            ${sample}.merged.krakenuniq.microbiome.fasta.human_synthetic.txt \\
            ${sample}.merged.krakenuniq.info_collection.flt.tmp \\
            > ${sample}.merged.krakenuniq.info_collection.flt
    else
        mv ${sample}.merged.krakenuniq.microbiome.fasta.tmp ${sample}.merged.krakenuniq.microbiome.fasta
        mv ${sample}.merged.krakenuniq.info_collection.flt.tmp ${sample}.merged.krakenuniq.info_collection.flt
    fi

    rm -f ${sample}.merged.krakenuniq.microbiome.fasta.human_synthetic.txt \\
          ${sample}.merged.krakenuniq.microbiome.fasta.tmp \\
          ${sample}.merged.krakenuniq.info_collection.flt.tmp
    """
}
