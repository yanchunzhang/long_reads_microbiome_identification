rule unmapped_analysis:
    input:
        "{sample}/{sample}.bam"
    output:
        "{sample}/{sample}.after_t2t.unmapped.fasta.gz"
    threads: 10
    resources:
        mem_mb=8000
    log:
        "logs/{sample}.unmapped.log"
    shell:
        """
        sh {config[scriptsdir]}/unmapped_analysis.sh \
        {wildcards.sample} {input} {wildcards.sample} {threads} \
        {config[scriptsdir]} {config[t2t_ref]} > {log} 2>&1
        """
