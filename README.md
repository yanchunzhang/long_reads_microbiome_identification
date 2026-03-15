Overview of the workflow.
This workflow identifies microbial reads from long-read sequencing data and evaluates their fragment length relative to host DNA. Starting from a BAM file with reads mapped to the human reference genome (hg38), unmapped reads are extracted and further filtered by alignment to the T2T reference genome. The remaining reads are classified using KrakenUniq with the MicrobialDB database, and candidate microbial reads are validated by megablast against the NCBI nt database. For each read, the best BLAST hit with alignment coverage greater than 0.5 is retained. Microbial reads are summarized at the genus level (≥5 reads per genus), and the Median(L)adj ratio is calculated as the ratio between the median microbial read length and the median human read length derived from samtools stats.

1. Human read filtering
Starting from a BAM file with reads mapped to the hg38 human reference genome, unmapped reads are extracted using samtools. These reads are then aligned to the T2T reference genome using minimap2, and reads that remain unmapped after both steps are retained as candidate non-human reads.
2. Human read length distribution
samtools stats is run on the original BAM file to obtain the read length distribution of human-mapped reads, which is used to calculate the median human read length.
3. Preliminary microbial classification
The candidate non-human reads are classified using KrakenUniq with the MicrobialDB database (downloaded from https://benlangmead.github.io/aws-indexes/k2) to obtain preliminary microbial taxonomic assignments.
4. Microbial read extraction and validation
Reads classified as microbial by KrakenUniq are extracted and aligned to the NCBI nt database using megablast. For each read, the best BLAST hit with the largest alignment coverage is selected. Only hits with alignment coverage > 0.5 are retained.
5. Genus-level fragment length analysis
Microbial reads are grouped at the genus level, considering only genera with at least five reads. For each genus, the median microbial read length is calculated.
6. Median(L)adj calculation
The Median(L)adj ratio is computed as:

**Median(L)adj** = median microbial read length (per genus) / median human read length
This metric provides a normalized measure of microbial fragment length relative to host DNA.


Reference
https://www.biorxiv.org/content/10.64898/2026.02.02.703393v1
