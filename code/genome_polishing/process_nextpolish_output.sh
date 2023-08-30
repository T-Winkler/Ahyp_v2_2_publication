#!/bin/bash

# This script takes the output of NextPolish and reproduces the complete genome file, in the same order as it was before (adding back ambiguous N bases).
# (see master_thesis/code/process_nextpolish_output.sh for more information about the input file preparation)

# Setup
NPOUT=data/NextPolish/output/
NPPROCESSED=data/NextPolish/processed/

mkdir -p "$NPPROCESSED"



# rename the nextpolish genome fasta
cut -f1,2 -d'_' "$NPOUT"genome.nextpolish.fa > "$NPPROCESSED"genome.nextpolish.renamed.fa

# index the prefiltered fasta file for use in R
samtools faidx data/NextPolish/input/out.prefiltered.renamed.txt

# filter the prefiltered file for everything that is not in Nextpolish file:
/home/tom/Documents/tools/bbmap/filterbyname.sh in=data/NextPolish/input/out.prefiltered.renamed.txt \
	names="$NPPROCESSED"genome.nextpolish.renamed.fa \
	out="$NPPROCESSED"prefilter_not_in_Nextpolish.fa

# use the helper Rscript to get the remaining headers without fasta sequences:
Rscript code/genome_polishing/helper_script.R

# Concatenate everything including the R output
cat "$NPPROCESSED"genome.nextpolish.renamed.fa \
	"$NPPROCESSED"prefilter_not_in_Nextpolish.fa \
	"$NPPROCESSED"header_without_sequence.fa > "$NPPROCESSED"combined.fa

# Linearize fasta file, this also adds a line after the empty fasta sequences in the end:
LC_ALL=C awk -v RS=">" -v FS="\n" -v ORS="\n" -v OFS="" '$0 {$1=">"$1"\n"; print}' "$NPPROCESSED"combined.fa > "$NPPROCESSED"combined.linear.fa

# rename header file by removing trailing underscore character of Contigs:
sed 's/quiver_/quiver/' data/NextPolish/input/out.headers.txt > "$NPPROCESSED"out.header.renamed.txt

# order file:
ORDER="$NPPROCESSED"out.header.renamed.txt
SORT="$NPPROCESSED"combined.linear.fa
OUT="$NPPROCESSED"combined.linear.sorted.fa

while read ID; do grep -w -A1 "$ID" $SORT; done < $ORDER > $OUT

# last step: remove "spl" lines and then all empty lines:
sed '/spl/d' "$NPPROCESSED"combined.linear.sorted.fa | awk 'NF' > "$NPPROCESSED"combined.linear.sorted.nosplit.fa

# normalize sequence length per line:
/home/tom/Documents/tools/gatk-4.2.5.0/gatk NormalizeFasta -I "$NPPROCESSED"combined.linear.sorted.nosplit.fa -O "$NPPROCESSED"combined.linear.sorted.nosplit.normalized.fa

# final output is data/NextPolish/processed/Ahypochondriacus_2.2_polished.fasta file
mv "$NPPROCESSED"combined.linear.sorted.nosplit.normalized.fa "$NPPROCESSED"Ahypochondriacus_2.2_polished.fasta

# remove intermediate output files:
rm "$NPPROCESSED"combined*
rm "$NPPROCESSED"genome*
rm "$NPPROCESSED"header*
rm "$NPPROCESSED"out*
rm "$NPPROCESSED"prefilter*

# copy to final output directory
mkdir -p polished_reference_genome/polished_genome_annotation/assembly/
cp "$NPPROCESSED"Ahypochondriacus_2.2_polished.fasta polished_reference_genome/polished_genome_annotation/assembly/
