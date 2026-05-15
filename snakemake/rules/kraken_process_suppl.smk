rule kraken_process_suppl:
    input:
        "{sample}.suppl.krakenuniq"
    output:
        fasta="{sample}.suppl.krakenuniq.microbiome.fasta",
        info="{sample}.suppl.krakenuniq.info_collection.flt"
    threads: 1
    resources:
        mem_mb=8000
    log:
        "logs/{sample}.kraken_process_suppl.log"
    shell:
        """
        sh {config[scriptsdir]}/post_kraken_filter.sh \
        {wildcards.sample}.suppl {config[kraken_db_suppl]} > {log} 2>&1
        """
