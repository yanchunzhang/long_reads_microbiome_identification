rule mapped_human_stats:
    input:
        bam="{sample}.bam"
    output:
        stats="{sample}.bam.mapped_human_reads_only.stats.txt"
    threads:
        config.get("samtools_stats_threads", 10)
    resources:
        mem_mb=config.get("samtools_stats_mem_mb", 6000)
    log:
        "logs/{sample}.mapped_human_stats.log"
    shell:
        r"""
        set -euo pipefail

        sh {config[scriptsdir]}/samtools_stats.mapped_human_reads_only.sh \
          {input.bam} \
          {threads} > {log} 2>&1
        """