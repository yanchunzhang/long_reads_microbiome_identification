rule krakenuniq_suppl:
    input:
        "{sample}.after_t2t.unmapped.fasta.gz"
    output:
        "{sample}.suppl.krakenuniq"
    threads: config.get("threads", {}).get("krakenuniq", 4)
    resources:
        mem_mb=config.get("resources", {}).get("krakenuniq_suppl_mem_mb", 60000)
    log:
        "logs/{sample}.kraken_suppl.log"
    shell:
        """
        sh {config[scriptsdir]}/krakenuniq.single.sh \
        {config[kraken_db_suppl]} \
        {input} \
        {wildcards.sample}.suppl \
        {threads} > {log} 2>&1
        """
