#!/bin/bash

# Script is used for the identification of betalain and flavonoid pathway genes in the A. hypochondriacus genome annotation version 2.2

# identification of betalains by blast
# create blast database of annotated protein sequences
# source of betalain pathway protein sequences is described in manuscript
mkdir -p data/annotation_analysis/betalains/DB
mkdir -p data/annotation_analysis/flavonoids

makeblastdb -in polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.prot.fasta \
	-out data/annotation_analysis/betalains/DB/blast_db \
	-logfile data/annotation_analysis/betalains/blast_db.log \
	-dbtype prot

# run protein blast of described betalain pathway genes against blast database
blastp -query data/annotation_analysis/betalains/pathway.fasta \
	-db data/annotation_analysis/betalains/DB/blast_db \
	-outfmt 7 \
	-out data/annotation_analysis/betalains/pathway_against_protein.out \
	-qcov_hsp_perc 80


# identification of flavonoids using KIPEs
# run KIPEs:
python /home/tom/Documents/tools/KIPEs/KIPEs3.py --baits /home/tom/Documents/tools/KIPEs/flavonoid_baits/ \
	--positions /home/tom/Documents/tools/KIPEs/flavonoid_residues/ \
	--out data/annotation_analysis/flavonoids/ \
	--subject polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.prot.fasta \
	--seqtype pep \
	--cpus 6
