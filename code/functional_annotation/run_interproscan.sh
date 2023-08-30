#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 168:00:00
#SBATCH -J interpro
#SBATCH -o logs/functional_annotation/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=32gb
#SBATCH --mail-type=ALL

# run on cheops1
# this script is used to run interproscan on the protein output of the manually corrected amaranth annotation
# installed on scratch due to the size of the installation
# installed interproscan using the following commands based on https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html
# downloaded with: wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.56-89.0/interproscan-5.56-89.0-64-bit.tar.gz
# checked md5sum after download
# extracted tarball and indexed hmm models before the first run using the following command: python3 initial_setup.py

# load required modules
module load openjdk/11.0.2

# set variables
# interpro directory with executable shell script
INTERPRODIR=/scratch/twinkle1/interproscan/interproscan-5.56-89.0/
# output directory for the annotation run
INTERPROOUT=data/functional_annotation/interproscan/
# input amino acid fasta file
INTERPROIN=polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.prot.fasta

# main
mkdir "$INTERPROOUT"
cd "$INTERPRODIR"

# dp:disables online lookup, f:output formats, iprlookup: for goterms and pa(thway) matching, 18 cpus as the main application also needs always 1 thread
./interproscan.sh -i "$INTERPROIN" \
	-dp \
	-f tsv,xml,gff3 \
	-d "$INTERPROOUT" \
	-iprlookup \
	-goterms \
	-pa \
	-cpu 18
