rule process_blast:
    input:
        blast="{sample}/{sample}.blast.txt"
    output:
        processed="{sample}/{sample}.blast.processed.txt",
        equivalent_hits="{sample}/{sample}.blast.equivalent_top_hits.tsv",
        equivalent_lca="{sample}/{sample}.blast.equivalent_top_hits.lca.tsv"
    threads: config.get("threads", {}).get("megablast", 10)
    resources:
        mem_mb=config.get("resources", {}).get("blast_process_mem_mb", 8000)
    params:
        taxonkit_data_dir=config.get("taxonkit_data_dir", "/sc/arion/projects/schzrnas/zhangy40/softwares/taxdump")
    log:
        "logs/{sample}.process_blast.log"
    shell:
        r"""
        set -o pipefail

        python {config[scriptsdir]}/blast_result_process.mt.py \
          --input {input.blast} \
          --output {output.processed} \
          --audit-output {output.equivalent_hits} \
          --threads {threads} \
          > {log} 2>&1

        (head -n 1 {output.equivalent_hits} | \
           awk 'BEGIN {{FS=OFS="\t"}} {{print $0,"representative_lca_taxid","representative_species_name","representative_species_taxid","equivalent_hit_lca_taxid","equivalent_hit_lca_lineage","equivalent_hit_lca_name","equivalent_hit_lca_rank","shared_species_name","shared_species_taxid","species_supported_by_equivalent_hits","representative_species_requires_tie_break"}}'; \
         tail -n +2 {output.equivalent_hits} | \
           taxonkit lca --data-dir {params.taxonkit_data_dir} -i 3 -s ';' -U -D | \
           taxonkit reformat --data-dir {params.taxonkit_data_dir} -I 10 -f '{{s}}' -t -r '' -R '' -T | \
           taxonkit lca --data-dir {params.taxonkit_data_dir} -i 7 -s ';' -U -D | \
           taxonkit lineage --data-dir {params.taxonkit_data_dir} -i 13 -n -r | \
           taxonkit reformat --data-dir {params.taxonkit_data_dir} -I 13 -f '{{s}}' -t -r '' -R '' -T | \
           awk 'BEGIN {{FS=OFS="\t"}} {{species_ok=($18!=""); tie_break=($5>1 && $12!="" && !species_ok); print $0,(species_ok ? "yes" : "no"),(tie_break ? "yes" : "no")}}') \
          > {output.equivalent_lca} 2>> {log}
        """

rule annotate_blast_lengths:
    input:
        info="{sample}/{sample}.merged.krakenuniq.info_collection.flt" if USE_SUPPL_DB else "{sample}/{sample}.krakenuniq.info_collection.flt",
        processed="{sample}/{sample}.blast.processed.txt"
    output:
        add_length="{sample}/{sample}.blast.processed.add_length.txt",
        other_microbiome="{sample}/{sample}.blast.microbiome.txt"
    resources:
        mem_mb=config.get("resources", {}).get("blast_process_mem_mb", 8000)
    log:
        "logs/{sample}.annotate_blast_lengths.log"
    shell:
        r"""
        set -o pipefail

        awk 'NR==FNR {{a[$3]=$8; next}} ($1 in a) {{print $0"\t"a[$1]}}' \
          {input.info} {input.processed} | \
        awk '{{print $0, $4/$5}}' | \
        sed 's/ /\t/g' | \
        taxonkit lca -i 3 -s ';' -U -D | \
        taxonkit reformat -I 7 -F -P | \
        awk 'BEGIN {{FS=OFS="\t"}} {{lineage=$8; $7=lineage; NF=7; print}}' | \
        sed 's/ /_/g' | \
        sort -k7,7 -k3,3 > {output.add_length} 2> {log}

        awk '$6>0.5 && !/k__unclass/ && !/g__unclass/ && /k__/ && \
             !/k__Metazoa/ && \
             (!/k__Euka/ || /(p__Ascomycota|p__Basidiomycota|p__Mucoromycota|p__Chytridiomycota)/)' \
          {output.add_length} > {output.other_microbiome} 2>> {log}
        """
