#!/bin/bash -l
#SBATCH -D /scratch/twinkle1/repeatmodeler/
#SBATCH -t 200:00:00
#SBATCH -J rmodeler
#SBATCH -o /home/twinkle1/master_thesis/logs/repeatmodeler/mappingLog-%j.txt
#SBATCH --error /home/twinkle1/master_thesis/logs/repeatmodeler/errorLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=42gb
#SBATCH --mail-user=twinkle1@smail.uni-koeln.de
#SBATCH --mail-type=ALL


# run on cheops1
# this script is used to run repeatmodeler on the newly polished reference assembly


# load modules
module load repeatmodeler/2.0.1

# create database directory
mkdir -p data/repeatmodeler/database/

### Main
# increased number of tasks, each rmblast job will take 4 threads
# for 20 threads, the pa setting while using rmblast should be 5

# Create database
BuildDatabase -name data/repeatmodeler/database/polished \
	data/NextPolish/processed/Ahypochondriacus_2.2_polished.fasta

# run Repeatmodeler
RepeatModeler -database data/repeatmodeler/database/polished \
	-pa 5 \
	-LTRStruct

# reclassify identified repeats
mkdir -p data/repeatmasking/reclassification/

# repeatmasker version 4.1.5
RepeatClassifier -consensi data/repeatmasking/reclassification/consensi.fa -stockholm data/repeatmasking/reclassification/families.stk
