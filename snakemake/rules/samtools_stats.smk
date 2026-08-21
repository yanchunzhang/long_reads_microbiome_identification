rule mapped_human_stats:
    input:
        bam="{sample}/{sample}.bam"
    output:
        stats="{sample}/{sample}.bam.mapped_human_reads_only.stats.txt.gz"
    threads:
        config.get("samtools_stats_threads", 10)
    resources:
        mem_mb=config.get("samtools_stats_mem_mb", 6000)
    log:
        "logs/{sample}.mapped_human_stats.log"
    shell:
        r"""
        set -euo pipefail
        source {config[scriptsdir]}/../../lib/hpc_modules.sh
        load_tool samtools samtools/1.21
        require_tools samtools

        samtools view -@ {threads} {input.bam} -h -F4 | \
          samtools stats 2> {log} | \
          gzip > {output.stats}
        """