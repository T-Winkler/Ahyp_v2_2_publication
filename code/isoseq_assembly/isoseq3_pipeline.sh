#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 12:00:00
#SBATCH -J isoseq3
#SBATCH -o logs/isoseq3/mappingLog-%j.txt
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --array 0-6


source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate isoseq3


echo "$SLURM_ARRAY_TASK_ID"

# following https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md from step 4 on using ccs.bam and primer.fa
# removal of primers and barcodes already done
# step 4, removal of polyA tail and artificial concatemers

PROCESSING=raw_data/isoseq_raw_reads/processed/
mkdir -p "PROCESSING"

SOURCE_DIR=raw_data/isoseq_raw_reads/reads
FILES=("$SOURCE_DIR"/*bam)

PRIMER_SOURCE=raw_data/isoseq_raw_reads/primers
PRIMER_FILES=("$PRIMER_SOURCE"/*fasta)

REFINEIN="$SOURCE_DIR"/"${FILES["${SLURM_ARRAY_TASK_ID}"]}"
PRIMERS="$PRIMER_SOURCE"/"${PRIMER_FILES["${SLURM_ARRAY_TASK_ID}"]}"
REFINEOUT="$PROCESSING""${FILES["${SLURM_ARRAY_TASK_ID}"]/.bam/flnc.bam}"

isoseq3 refine $REFINEIN $PRIMERS $REFINEOUT --require-polya

# flnc reads were created using the following command, removing artificial concatemers and filtering out all genes without Poly-A tail
# using respective primers of the different datasets
#isoseq3 refine $REFINEIN $PRIMERS $REFINEOUT --require-polya


# merge all flnc reads into a single file
samtools merge "$ISOSEQOUT"combined.merged_flnc.bam /home/twinkle1/isoseq_raw_reads/*bam


# clustering of identical reads

IN="data/isoseq3_pipeline/combined.merged_flnc.bam"
OUT="data/isoseq3_pipeline/combined.merged_clustered.bam"

isoseq3 cluster	"$IN" "$OUT" --verbose --use-qvs
