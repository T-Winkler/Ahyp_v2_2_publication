#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 160:00:00
#SBATCH -J rmodeler
#SBATCH -o /projects/ag-stetter/twinkle/projects/Ahyp_v2_2/logs/repeatmasker/mappingLog-%j.txt
#SBATCH --error /projects/ag-stetter/twinkle/projects/Ahyp_v2_2/logs/repeatmasker/errorLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=42gb
#SBATCH --mail-type=ALL


# run on cheops1
# this script is used to run repeatmasker usuing the generated repeatmodeler repeatdatabase

# load modules
module load repeatmasker/4.1.1

# create database directory
RMOUT=data/repeatmasking/repeatmasker

mkdir -p $RMOUT

### Main
# watch out with the processor setting, each rmblast job will take 4 threads
# for 20 threads, the pa setting while using rmblast should be 5

# Converted the polished reference genome fasta to all capital letters as preparation to repeatmasking:
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' data/NextPolish/processed/Ahypochondriacus/V2_2/Ahypochondriacus_2.2_polished.softmasked.fasta \
	> "$RMOUT"/Ahypochondriacus_2.2_polished.capital.fasta

### run RepeatMasker
# lib = repeat database created with repeatmodeler
RepeatMasker -lib data/repeatmasking/repeatmodeler/consensi.fa.classified \
	-pa 5 \
	-small \
	-e rmblast \
	-gff \
	-dir "$RMOUT" \
	"$RMOUT"/Ahypochondriacus_2.2_polished.capital.fasta

mkdir "$RMOUT"/output

# convert to softmasked fasta
module load bedtools/2.29.2
bedtools maskfasta \
	-fi "$RMOUT"/Ahypochondriacus_2.2_polished.capital.fasta \
	-bed "$RMOUT"/Ahypochondriacus_2.2_polished.capital.fasta.out.gff \
	-soft \
	-fo "$RMOUT"/output/Ahypochondriacus_2.2_polished.softmasked.fasta

# more detailed summary:
buildSummary.pl data/repeatmasking/repeatmasker/Ahypochondriacus_2.2_polished.capital.fasta.out > data/repeatmasking/repeatmasker/Ahypochondriacus_2.2_polished.capital.fasta.buildSummary
