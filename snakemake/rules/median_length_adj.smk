rule median_length_adj:
    input:
        stats="{sample}/{sample}.bam.mapped_human_reads_only.stats.txt.gz",
        microbe="{sample}/{sample}.blast.microbiome.txt"

    output:
        median_L_adj="{sample}/{sample}.median_l_adj.txt",
        sum="{sample}/{sample}.microbiome.sum_by_length_per_genus.txt",
        human_median="{sample}/{sample}.human_median_length.tsv"

    log:
        "logs/{sample}.median_length.log"

    shell:
        """
        python {config[scriptsdir]}/median_length_adj.py \
            --stats {input.stats} \
            --microbe {input.microbe} \
            --sum_out {output.sum} \
            --gt5_out {output.median_L_adj} \
            --human_median_out {output.human_median} \
            --sample {wildcards.sample} > {log} 2>&1
        """
