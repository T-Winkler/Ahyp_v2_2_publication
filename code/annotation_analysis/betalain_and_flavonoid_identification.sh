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

# create blast db
mkdir data/annotation_analysis/flavonoids/blast_db
makeblastdb -in polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta -out data/annotation_analysis/flavonoids/blast_db/db -dbtype nucl

# blast searches for unidentified candidate genes using KIPEs bait sequences
tblastn -query data/annotation_analysis/flavonoids/KIPEs/blast_query/ANS.fasta -db data/annotation_analysis/flavonoids/blast_db/db -outfmt 7 > data/annotation_analysis/flavonoids/ANS_blast.out
tblastn -query data/annotation_analysis/flavonoids/KIPEs/blast_query/ANR.fasta -db data/annotation_analysis/flavonoids/blast_db/db -outfmt 7 > data/annotation_analysis/flavonoids/ANR_blast.out
tblastn -query data/annotation_analysis/flavonoids/KIPEs/blast_query/F3-5-H.fasta -db data/annotation_analysis/flavonoids/blast_db/db -outfmt 7 > data/annotation_analysis/flavonoids/F3-5-H_blast.out
tblastn -query data/annotation_analysis/flavonoids/KIPEs/blast_query/FNS1.fasta -db data/annotation_analysis/flavonoids/blast_db/db -outfmt 7 > data/annotation_analysis/flavonoids/FNS1_blast.out
tblastn -query data/annotation_analysis/flavonoids/KIPEs/blast_query/CHI2.fasta -db data/annotation_analysis/flavonoids/blast_db/db -outfmt 7 > data/annotation_analysis/flavonoids/CHI2_blast.out

# exonerate protein alignments
# prepare fasta file of KIPEs bait sequences
sed 's/%_//' data/annotation_analysis/flavonoids/KIPEs/blast_query/ANS.fasta > data/annotation_analysis/flavonoids/ANS_fixed.fasta
# run exonerate
exonerate --model protein2genome --percent 55 --showvulgar no --showalignment no --showtargetgff yes --query data/annotation_analysis/flavonoids/ANS_fixed.fasta --target polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta > data/annotation_analysis/flavonoids/ANS_exonerate.gff

sed 's/%_//' data/annotation_analysis/flavonoids/KIPEs/blast_query/ANR.fasta > data/annotation_analysis/flavonoids/ANR_fixed.fasta
exonerate --model protein2genome --percent 55 --showvulgar no --showalignment no --showtargetgff yes --query data/annotation_analysis/flavonoids/ANR_fixed.fasta --target polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta > data/annotation_analysis/flavonoids/ANR_exonerate.gff

sed 's/%_//' data/annotation_analysis/flavonoids/KIPEs/blast_query/F3-5-H.fasta > data/annotation_analysis/flavonoids/F3-5-H_fixed.fasta
exonerate --model protein2genome --percent 55 --showvulgar no --showalignment no --showtargetgff yes --query data/annotation_analysis/flavonoids/F3-5-H_fixed.fasta --target polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta > data/annotation_analysis/flavonoids/F3-5-H_exonerate.gff

sed 's/%_//' data/annotation_analysis/flavonoids/KIPEs/blast_query/FNS1.fasta > data/annotation_analysis/flavonoids/FNS1_fixed.fasta
exonerate --model protein2genome --percent 55 --showvulgar no --showalignment no --showtargetgff yes --query data/annotation_analysis/flavonoids/FNS1_fixed.fasta --target polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta > data/annotation_analysis/flavonoids/FNS1_exonerate.gff

sed 's/%_//' data/annotation_analysis/flavonoids/KIPEs/blast_query/CHI2.fasta > data/annotation_analysis/flavonoids/CHI2_fixed.fasta
exonerate --model protein2genome --percent 55 --showvulgar no --showalignment no --showtargetgff yes --query data/annotation_analysis/flavonoids/CHI2_fixed.fasta --target polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta > data/annotation_analysis/flavonoids/CHI2_exonerate.gff
