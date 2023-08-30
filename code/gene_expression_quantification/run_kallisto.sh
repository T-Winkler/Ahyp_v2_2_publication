#!/bin/bash -l

# Beforehand:
#Short read data downloaded from SRA using the following accession/run numbers:

#SRA Accession numbers:
#Floral tissue: SRX722058 SRR1598911
#Leaf tissue: SRX722059 SRR1598912
#Root tissue: SRX722060 SRR1598913
#Stem tissue: SRX722057 SRR1598910
#Water stressed tissue sample: SRX722061 SRR1598914
#Immature seeds: SRX722056 SRR1598909
#Mature seeds: SRX722063 SRR1598916
#Green Cotyledone: SRX722062 SRR1598915


# index the transcriptome
KALINDEX=data/gene_expression_quantification/kallisto_quant/index
mkdir -p $KALINDEX

/home/tom/Documents/tools/kallisto/kallisto index -i "$KALINDEX"/index polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.cds.fasta


##### Perform quantification
# create array of read fastq files (R1 only):
SOURCE_DIR=raw_data/Clouse_short_reads/
FILES=("$SOURCE_DIR"SRR*_1.fastq.gz)
OUTDIR=data/gene_expression_quantification/kallisto_quant/
TISSUE_NAMES=("young_seed" "stem" "flower" "leaf" "root" "water_stressed" "cotyledones" "mature_seed")

# kallisto after indexing
for (( i=0; i<=7; i++))
do
	/home/tom/Documents/tools/kallisto/kallisto quant -i "$KALINDEX"/index -o "$OUTDIR""${TISSUE_NAMES[$i]}" --bias --plaintext -t 6 --verbose "${FILES[$i]}" "${FILES[$i]/_1.fastq.gz/_2.fastq.gz}"
done
