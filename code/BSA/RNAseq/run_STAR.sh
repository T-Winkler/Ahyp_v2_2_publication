#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 1:00:00
#SBATCH -J STAR
#SBATCH -o /home/twinkle1/projects/Ahyp_v2_2/logs/flower_color_mapping/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=8gb
#SBATCH --array 0-3


module load star/2.7.8a

# create array of read fastq files (R1 only):
SOURCE_DIR=raw_data/BSA_rnaseq/
FILES=("$SOURCE_DIR"/*_1P.fq.gz)
OUTPUTDIR=data/BSA/RNAseq/STAR_flower_mappings

# change directory of outprefix
OUTPREFIX1=("${FILES["${SLURM_ARRAY_TASK_ID}"]/$SOURCE_DIR/$OUTPUTDIR}")
# change suffix of outprefix by removing the file extension etc.
OUTPREFIX="${OUTPREFIX1/trimmed_1P.fq.gz/}"


# run STAR after genome index creation
mkdir -p "$OUTPUTDIR"

STAR --runThreadN 8 \
	--runMode alignReads \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir /scratch/twinkle1/STAR_flower_index \
	--outFileNamePrefix "$OUTPREFIX" \
	--readFilesCommand zcat \
	--readFilesIn "${FILES["${SLURM_ARRAY_TASK_ID}"]}" "${FILES["${SLURM_ARRAY_TASK_ID}"]/_1P.fq.gz/_2P.fq.gz}"
