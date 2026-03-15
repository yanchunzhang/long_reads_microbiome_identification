rule kraken_process:
    input:
        "{sample}.krakenuniq"
    output:
        fasta="{sample}.krakenuniq.microbiome.fasta",
	info="{sample}.krakenuniq.info_collection.flt"
    threads: 1
    resources:
        mem_mb=8000
    log:
        "logs/{sample}.kraken_process.log"
    shell:
        """
        sh {config[scriptsdir]}/post_kraken_filter.sh \
        {wildcards.sample} > {log} 2>&1
        """
