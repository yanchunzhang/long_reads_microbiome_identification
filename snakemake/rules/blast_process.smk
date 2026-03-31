rule process_blast:
    input:
        blast="{sample}.blast.txt"
    output:
        processed="{sample}.blast.processed.txt"
    threads: config.get("threads", {}).get("megablast", 10)
    resources:
        mem_mb=config.get("resources", {}).get("blast_process_mem_mb", 8000)
    log:
        "logs/{sample}.process_blast.log"
    shell:
        r"""
        python {config[scriptsdir]}/blast_result_process.mt.py \
          --input {input.blast} \
          --output {output.processed} \
          --threads {threads} \
          > {log} 2>&1
        """

rule annotate_blast_lengths:
    input:
        info="{sample}.krakenuniq.info_collection.flt",
        processed="{sample}.blast.processed.txt"
    output:
        add_length="{sample}.blast.processed.add_length.txt",
        other_microbiome="{sample}.blast.microbiome.txt"
    resources:
        mem_mb=config.get("resources", {}).get("blast_process_mem_mb", 8000)
    log:
        "logs/{sample}.annotate_blast_lengths.log"
    shell:
        r"""
        set -o pipefail

        awk 'NR==FNR {{a[$3]=$8; next}} ($1 in a) {{print $0"\t"a[$1]}}' \
          {input.info} {input.processed} | \
        awk '{{print $0, $4/$5}}' | \
        sed 's/ /\t/g' | \
        taxonkit reformat -I 3 -F -P | \
        sed 's/ /_/g' | \
        sort -k7,7 -k3,3 > {output.add_length} 2> {log}

        awk '$6>0.5 && !/k__unclass/ && !/g__unclass/ && !/k__Euka/ && /k__/' \
          {output.add_length} > {output.other_microbiome} 2>> {log}
        """
