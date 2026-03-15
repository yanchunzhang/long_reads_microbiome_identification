rule median_length_adj:
    input:
        stats="{sample}.bam.mapped_human_reads_only.stats.txt",
        microbe="{sample}.blast.microbiome.txt"

    output:
        median_L_adj="{sample}.median_l_adj.txt",
        sum="{sample}.microbiome.sum_by_length_per_genus.txt"

    log:
        "logs/{sample}.median_length.log"

    shell:
        """
        python {config[scriptsdir]}/median_length_adj.py \
            --stats {input.stats} \
            --microbe {input.microbe} \
            --sum_out {output.sum} \
            --gt5_out {output.median_L_adj} \
            --sample {wildcards.sample} > {log} 2>&1
        """
