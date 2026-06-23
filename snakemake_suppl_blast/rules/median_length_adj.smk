rule median_length_adj_virus_fungi:
    input:
        stats="{sample}.bam.mapped_human_reads_only.stats.txt.gz",
        microbe="{sample}/{sample}.blast.virus_fungi.txt"

    output:
        median_L_adj="{sample}/{sample}.virus_fungi.median_l_adj.txt",
        sum="{sample}/{sample}.virus_fungi.sum_by_length_per_genus.txt"

    log:
        "logs/{sample}.median_length_virus_fungi.log"

    shell:
        """
        python {config[scriptsdir]}/median_length_adj.py \
            --stats {input.stats} \
            --microbe {input.microbe} \
            --sum_out {output.sum} \
            --gt5_out {output.median_L_adj} \
            --sample {wildcards.sample} > {log} 2>&1
        """


rule median_length_adj:
    input:
        stats="{sample}.bam.mapped_human_reads_only.stats.txt.gz",
        microbe="{sample}/{sample}.blast.microbiome.txt"

    output:
        median_L_adj="{sample}/{sample}.median_l_adj.txt",
        sum="{sample}/{sample}.microbiome.sum_by_length_per_genus.txt"

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
