// modules/blast_process.nf
// Two post-BLAST processing steps:
//   PROCESS_BLAST          – multi-threaded BLAST hit processing (Python)
//   ANNOTATE_BLAST_LENGTHS – join KrakenUniq read-length info, reformat with
//                            taxonkit, and filter to microbial kingdom hits
//   FILTER_ONT_ARTIFACTS   – merge ONT construct intervals and remove reads
//                            with technical coverage at or above the threshold
//
// Mirrors Snakemake rules: process_blast, annotate_blast_lengths


// ── 1. Process raw BLAST output ───────────────────────────────────────────
process PROCESS_BLAST {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(blast)

    output:
    tuple val(sample), path("${sample}.blast.processed.txt"), emit: processed
    tuple val(sample), path("${sample}.blast.equivalent_top_hits.tsv"), emit: equivalent_hits
    tuple val(sample), path("${sample}.blast.equivalent_top_hits.lca.tsv"), emit: equivalent_lca

    script:
    """
    set -o pipefail

    python ${params.scriptsdir}/blast_result_process.mt.py \\
        --input   ${blast} \\
        --output  ${sample}.blast.processed.txt \\
        --audit-output ${sample}.blast.equivalent_top_hits.tsv \\
        --threads ${task.cpus}

    (head -n 1 ${sample}.blast.equivalent_top_hits.tsv | \\
       awk 'BEGIN {FS=OFS="\t"} {print \$0,"representative_lca_taxid","representative_species_name","representative_species_taxid","equivalent_hit_lca_taxid","equivalent_hit_lca_lineage","equivalent_hit_lca_name","equivalent_hit_lca_rank","shared_species_name","shared_species_taxid","species_supported_by_equivalent_hits","representative_species_requires_tie_break"}'; \\
     tail -n +2 ${sample}.blast.equivalent_top_hits.tsv | \\
       taxonkit lca --data-dir ${params.taxonkit_data_dir} -i 3 -s ';' -U -D | \\
       taxonkit reformat --data-dir ${params.taxonkit_data_dir} -I 10 -f '{s}' -t -r '' -R '' -T | \\
       taxonkit lca --data-dir ${params.taxonkit_data_dir} -i 7 -s ';' -U -D | \\
       taxonkit lineage --data-dir ${params.taxonkit_data_dir} -i 13 -n -r | \\
       taxonkit reformat --data-dir ${params.taxonkit_data_dir} -I 13 -f '{s}' -t -r '' -R '' -T | \\
       awk 'BEGIN {FS=OFS="\t"} {species_ok=(\$18!=""); tie_break=(\$5>1 && \$12!="" && !species_ok); print \$0,(species_ok ? "yes" : "no"),(tie_break ? "yes" : "no")}') \\
      > ${sample}.blast.equivalent_top_hits.lca.tsv
    """
}


// ── 2. Annotate hits with read-length info and filter to microbiome ───────
// Input channel carries: [sample, info_flt, processed_txt]
// (joined from KRAKEN_PROCESS.out.info + PROCESS_BLAST.out in main.nf)
//
// Emits two named output channels:
//   .add_length  → [sample, *.blast.processed.add_length.txt]
//   .microbiome_pre_filter → preliminary microbial calls; FILTER_ONT_ARTIFACTS
//                             turns these into *.blast.microbiome.txt

process ANNOTATE_BLAST_LENGTHS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(info), path(processed)

    output:
    tuple val(sample), path("${sample}.blast.processed.add_length.txt"), emit: add_length
    tuple val(sample), path("${sample}.blast.microbiome.pre_ont_filter.txt"), emit: microbiome_pre_filter

    script:
    // Note: \$ is required inside Nextflow """ blocks to pass a literal $
    // to the shell (Groovy interpolates ${...} but leaves \$ as $).
    """
    set -o pipefail

    awk 'NR==FNR {a[\$3]=\$8; next} (\$1 in a) {print \$0"\\t"a[\$1]}' \\
        ${info} ${processed} | \\
    awk '{print \$0, \$4/\$5}' | \\
    sed 's/ /\\t/g' | \\
    taxonkit lca -i 3 -s ';' -U -D | \\
    taxonkit reformat -I 7 -F -P | \\
    awk 'BEGIN {FS=OFS="\\t"} {lineage=\$8; \$7=lineage; NF=7; print}' | \\
    sed 's/ /_/g' | \\
    sort -k7,7 -k3,3 > ${sample}.blast.processed.add_length.txt

    awk '\$6>${params.blast_min_query_coverage} && !/k__unclass/ && !/g__unclass/ && /k__/ && \\
         !/k__Metazoa/ && \\
         (!/k__Euka/ || /(p__Ascomycota|p__Basidiomycota|p__Mucoromycota|p__Chytridiomycota)/)' \\
        ${sample}.blast.processed.add_length.txt > ${sample}.blast.microbiome.pre_ont_filter.txt
    """
}


// ── 3. Remove ONT construct-dominated reads ────────────────────────────────
// Qualifying barcode/adapter-hit intervals are merged, and a read is removed
// when their union reaches the configured fraction of the full read.
process FILTER_ONT_ARTIFACTS {
    tag "${sample}"

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(microbiome_pre_filter), path(query_fasta), path(adapter_fasta)

    output:
    tuple val(sample), path("${sample}.blast.microbiome.txt"), emit: microbiome
    tuple val(sample), path("${sample}.blast.ont_adapter_filter.audit.tsv"), emit: audit
    tuple val(sample), path("${sample}.blast.ont_adapter_filtered_out.txt"), emit: filtered_out
    tuple val(sample), path("${sample}.blast.ont_adapter_hits.tsv"), emit: adapter_hits

    script:
    def filterEnabled = params.filter_ont_adapters.toString().toLowerCase() != 'false'
    def disabled = filterEnabled ? '' : '--disabled'
    def enabledText = filterEnabled ? 'true' : 'false'
    """
    set -euo pipefail

    source ${projectDir}/../lib/hpc_modules.sh
    load_tool blastn blast/2.13.0+
    require_tools blastn

    if [[ "${enabledText}" == "false" ]]; then
        : > ${sample}.blast.ont_adapter_hits.tsv
    else
        blastn -task blastn-short \
          -query ${query_fasta} \
          -subject ${adapter_fasta} \
          -strand both \
          -word_size 7 \
          -perc_identity ${params.ont_adapter_min_identity} \
          -max_hsps 10 \
          -dust no \
          -soft_masking false \
          -evalue 1000 \
          -num_threads ${task.cpus} \
          -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
          -out ${sample}.blast.ont_adapter_hits.tsv
    fi

    python ${params.scriptsdir}/filter_ont_artifacts_after_blast.py \
      --microbiome ${microbiome_pre_filter} \
      --ont-hits ${sample}.blast.ont_adapter_hits.tsv \
      --output ${sample}.blast.microbiome.txt \
      --audit-output ${sample}.blast.ont_adapter_filter.audit.tsv \
      --filtered-output ${sample}.blast.ont_adapter_filtered_out.txt \
      --min-overlap ${params.ont_adapter_min_overlap} \
      --min-identity ${params.ont_adapter_min_identity} \
      --min-technical-fraction ${params.ont_adapter_min_technical_fraction} \
      ${disabled}
    """
}
