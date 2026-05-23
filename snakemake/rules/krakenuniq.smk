rule krakenuniq:
    input:
        "{sample}/{sample}.after_t2t.unmapped.fasta.gz"
    output:
        "{sample}/{sample}.krakenuniq"
    threads: 4
    resources:
        mem_mb=150000
    log:
        "logs/{sample}.kraken.log"
    shell:
        """
        sh {config[scriptsdir]}/krakenuniq.single.sh \
        {config[kraken_db]} \
        {input} \
        {wildcards.sample}/{wildcards.sample} \
        {threads} > {log} 2>&1
        """
