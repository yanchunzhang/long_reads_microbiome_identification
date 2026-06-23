rule kraken_process_suppl:
    input:
        "{sample}/{sample}.suppl.krakenuniq"
    output:
        fasta="{sample}/{sample}.suppl.krakenuniq.microbiome.fasta",
        info="{sample}/{sample}.suppl.krakenuniq.info_collection.flt"
    threads: 1
    resources:
        mem_mb=8000
    log:
        "logs/{sample}.kraken_process_suppl.log"
    shell:
        """
        sh {config[scriptsdir]}/post_kraken_filter.sh \
        {wildcards.sample}/{wildcards.sample}.suppl {config[kraken_db_suppl]} > {log} 2>&1
        """
