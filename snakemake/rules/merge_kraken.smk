rule merge_kraken:
    input:
        primary_fasta="{sample}.krakenuniq.microbiome.fasta",
        suppl_fasta="{sample}.suppl.krakenuniq.microbiome.fasta",
        primary_info="{sample}.krakenuniq.info_collection.flt",
        suppl_info="{sample}.suppl.krakenuniq.info_collection.flt"
    output:
        fasta="{sample}.merged.krakenuniq.microbiome.fasta",
        info="{sample}.merged.krakenuniq.info_collection.flt"
    log:
        "logs/{sample}.merge_kraken.log"
    shell:
        r"""
        set -euo pipefail

        # Merge FASTA — skip any read ID already seen (primary DB has priority)
        awk '/^>/ {{ if (seen[$0]) {{ skip=1 }} else {{ seen[$0]=1; skip=0; print }} }} \
             !/^>/ && !skip {{ print }}' \
            {input.primary_fasta} {input.suppl_fasta} > {output.fasta} 2>> {log}

        # Merge info collection — first occurrence wins (primary DB first)
        awk '!seen[$1]++' {input.primary_info} {input.suppl_info} > {output.info} 2>> {log}
        """
