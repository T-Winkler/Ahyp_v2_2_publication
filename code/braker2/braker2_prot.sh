#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 100:00:00
#SBATCH -J braker
#SBATCH -o logs/braker/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=84gb
#SBATCH --mail-user=twinkle1@smail.uni-koeln.de
#SBATCH --mail-type=ALL

#Braker2 protein input using the plant sequences from orthoDB (downloaded from: https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz) as well as the protein sequences from amaranthus cruentus (removed asterisks and space in fasta header)
#The downloaded dataset contains sequences from 117 embryophyte species.

#wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
#tar -xvf odb10_plants_fasta.tar.gz
# write all into the same file:
#cat plants/Rawdata/* > protein_db.fasta

# also add the Cruentus sequences:
#cat /projects/ag-stetter/twinkle/Amaranthus_cruentus/Amacr_pep_nospace.fa >> protein_db.fasta
# removed asterisk
#sed 's/\*//' protein_db.fasta > protein_db.fa

#-> in total 3536219 plant protein sequences
#Since I added the Cruentus sequences, the total number of species is now 118.


# run on cheops1
# this script is used to run braker2 in both RNAseq and protein mode

source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate /opt/rrzk/software/conda-envs/braker2

module load samtools/1.13

mkdir -p data/braker2/polished_prot

# run braker in RNAseq and protein mode:

braker.pl --AUGUSTUS_CONFIG_PATH=/home/twinkle1/tools/config/ \
	--epmode \
	--prot_seq=data/braker2/input/protein_db.fasta \
	--genome=polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta \
	--softmasking \
	--species=polished_prot \
	--cores=8 \
	--workingdir=data/braker2/polished_prot
