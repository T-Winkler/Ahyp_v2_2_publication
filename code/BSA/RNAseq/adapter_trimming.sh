#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 4:00:00
#SBATCH -J trimmomatic
#SBATCH -o /home/twinkle1/projects/Ahyp_v2_2/logs/flower_color_mapping/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=8gb
#SBATCH --array 0-3
#SBATCH --mail-user=twinkle1@smail.uni-koeln.de
#SBATCH --mail-type=ALL

# trim bulk segregant RNAseq data using trimmomatic with the specified adapter sequences

module load trimmomatic/0.39

# there are a total of 4 samples, array goes from 0-3

### MAIN

# run this part as array job
# create array of read fastq files (R1 only):
SOURCE_DIR=raw_data/BSA_rnaseq/
FILES=("$SOURCE_DIR"/*R1.fastq.gz)

# run trimmomatic, use 6 threads, taking advantage of the baseout function to name output files
# use the sequencing adapters send by the sequencing center as custom fasta file
# forward and reverse read adapters are indicated in the fasta file by the /1 and /2 suffixes

java -jar $TRIMMOMATIC/trimmomatic.jar PE \
	-threads 6 \
	"${FILES["${SLURM_ARRAY_TASK_ID}"]}" \
	"${FILES["${SLURM_ARRAY_TASK_ID}"]/R1.fastq.gz/R2.fastq.gz}" \
	-baseout "${FILES["${SLURM_ARRAY_TASK_ID}"]/R1.fastq.gz/trimmed.fq.gz}" \
	ILLUMINACLIP:raw_data/BSA_rnaseq/adapters/custom_adapters.fa:2:30:10
