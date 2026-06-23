#!/bin/bash
#BSUB -J create_test_bam
#BSUB -P acc_schzrnas
#BSUB -q express
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=8000]
#BSUB -W 2:00
#BSUB -o logs/create_test_bam.out
#BSUB -e logs/create_test_bam.err

module load samtools/1.21

BAM=/sc/arion/projects/schzrnas/zhangy40/intratumor_bacteria/snakemake/CO4_T1.bam
MICROBIOME=/sc/arion/projects/schzrnas/zhangy40/intratumor_bacteria/snakemake/CO4_T1.blast.microbiome.txt
OUTDIR=/sc/arion/projects/schzrnas/zhangy40/scripts/long_reads_microbiome_identification/long_reads_microbiome_identification/test
THREADS=4

mkdir -p $OUTDIR/tmp

# Step 1: Extract microbiome read IDs (strip leading @)
awk '{print substr($1,2)}' $MICROBIOME > $OUTDIR/tmp/microbiome_readids.txt
echo "Microbiome reads: $(wc -l < $OUTDIR/tmp/microbiome_readids.txt)"

# Step 2: Extract microbiome reads from original BAM
samtools view -@ $THREADS -N $OUTDIR/tmp/microbiome_readids.txt -b $BAM \
    -o $OUTDIR/tmp/microbiome.bam

# Step 3: Subset ~1% of mapped human reads (exclude unmapped -F4)
samtools view -@ $THREADS -F4 -s 42.01 -b $BAM \
    -o $OUTDIR/tmp/human_subset.bam

# Step 4: Merge and sort
samtools merge -@ $THREADS -f $OUTDIR/tmp/merged.bam \
    $OUTDIR/tmp/human_subset.bam \
    $OUTDIR/tmp/microbiome.bam

samtools sort -@ $THREADS $OUTDIR/tmp/merged.bam -o $OUTDIR/CO4_T1.test.bam
samtools index $OUTDIR/CO4_T1.test.bam

# Step 5: Report
echo "Test BAM size: $(ls -lh $OUTDIR/CO4_T1.test.bam | awk '{print $5}')"
echo "Test BAM reads: $(samtools view -c $OUTDIR/CO4_T1.test.bam)"

# Cleanup
rm -rf $OUTDIR/tmp
