import os

BLAST_SPLIT_NSEQ = int(config.get("blast_split_nseq", 5000))


checkpoint split_blast_query_fasta:
    input:
        fasta="{sample}.krakenuniq.microbiome.fasta"
    output:
        split_dir=directory("split_fasta/{sample}")
    params:
        nseq=BLAST_SPLIT_NSEQ
    log:
        "logs/{sample}.split_blast_query.log"
    shell:
        r"""
        set -euo pipefail

        rm -rf {output.split_dir}
        mkdir -p {output.split_dir}

        seqkit split2 \
          -s {params.nseq} \
          -O {output.split_dir} \
          {input.fasta} > {log} 2>&1

        i=0
        found=0
        for f in {output.split_dir}/*.fasta {output.split_dir}/*.fa; do
            [ -e "$f" ] || continue
            found=1
            new=$(printf "{output.split_dir}/{wildcards.sample}.part_%03d.fa" "$i")
            mv "$f" "$new"
            i=$((i+1))
        done

        if [ "$found" -eq 0 ]; then
            echo "ERROR: no split FASTA files were produced for {wildcards.sample}" >&2
            exit 1
        fi
        """
        

def get_split_fasta_files(wildcards):
    ckpt = checkpoints.split_blast_query_fasta.get(sample=wildcards.sample)
    split_dir = ckpt.output.split_dir

    files = sorted([
        os.path.join(split_dir, f)
        for f in os.listdir(split_dir)
        if f.endswith(".fa")
    ])

    if not files:
        raise ValueError(f"No split FASTA files found in {split_dir}")

    return files


def get_blast_chunk_outputs(wildcards):
    fasta_files = get_split_fasta_files(wildcards)
    outputs = []
    for f in fasta_files:
        chunk = os.path.basename(f).rsplit(".fa", 1)[0]
        outputs.append(f"blast_chunks/{wildcards.sample}/{chunk}.blast.txt")
    return outputs


rule megablast_chunk:
    input:
        fasta="split_fasta/{sample}/{chunk}.fa"
    output:
        blast=temp("blast_chunks/{sample}/{chunk}.blast.txt")
    threads:
        config.get("threads", {}).get("megablast", 8)
    resources:
        mem_mb=config.get("resources", {}).get("megablast_mem_mb", 30000)
    log:
        "logs/{sample}.{chunk}.megablast.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p blast_chunks/{wildcards.sample}

        sh {config[scriptsdir]}/megablast.sh \
          {input.fasta} \
          {config[blastdb]} \
          {output.blast} \
          {threads} blast > {log} 2>&1
        """


rule merge_blast_chunks:
    input:
        chunks=get_blast_chunk_outputs
    output:
        blast="{sample}.blast.txt"
    log:
        "logs/{sample}.merge_blast_chunks.log"
    run:
        with open(output.blast, "w") as out:
            for fn in sorted(input.chunks):
                with open(fn) as fh:
                    for line in fh:
                        out.write(line)