rule merge_kraken:
    input:
        primary_fasta="{sample}/{sample}.krakenuniq.microbiome.fasta",
        suppl_fasta="{sample}/{sample}.suppl.krakenuniq.microbiome.fasta",
        primary_info="{sample}/{sample}.krakenuniq.info_collection.flt",
        suppl_info="{sample}/{sample}.suppl.krakenuniq.info_collection.flt",
        primary_kraken="{sample}/{sample}.krakenuniq"
    output:
        fasta="{sample}/{sample}.merged.krakenuniq.microbiome.fasta",
        info="{sample}/{sample}.merged.krakenuniq.info_collection.flt"
    log:
        "logs/{sample}.merge_kraken.log"
    shell:
        r"""
        set -euo pipefail

        # Merge FASTA — skip any read ID already seen (primary DB has priority)
        awk '/^>/ {{ if (seen[$0]) {{ skip=1 }} else {{ seen[$0]=1; skip=0; print }} }} \
             !/^>/ && !skip {{ print }}' \
            {input.primary_fasta} {input.suppl_fasta} > {output.fasta}.tmp 2>> {log}

        # Merge info collection — first occurrence wins (primary DB first)
        awk '!seen[$1]++' {input.primary_info} {input.suppl_info} > {output.info}.tmp 2>> {log}

        # Remove reads the primary DB classified as human (9606) or synthetic (32630).
        # The suppl DB has no human sequences so it misclassifies human reads as
        # fungi/viruses; the primary DB correctly flags them. Pre-filtering here
        # prevents false positives and avoids blasting millions of human reads.
        awk '$3==9606 || $3==32630 {{print $2}}' {input.primary_kraken} \
            > {output.fasta}.human_synthetic.txt 2>> {log}

        if [[ -s "{output.fasta}.human_synthetic.txt" ]]; then
            seqkit grep -v -f {output.fasta}.human_synthetic.txt \
                {output.fasta}.tmp > {output.fasta} 2>> {log}
            awk 'NR==FNR {{seen[">"$1]=1; next}} !($1 in seen)' \
                {output.fasta}.human_synthetic.txt \
                {output.info}.tmp > {output.info} 2>> {log}
        else
            mv {output.fasta}.tmp {output.fasta}
            mv {output.info}.tmp {output.info}
        fi
        rm -f {output.fasta}.human_synthetic.txt {output.fasta}.tmp {output.info}.tmp
        """
