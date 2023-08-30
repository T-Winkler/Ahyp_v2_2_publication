#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 10:00:00
#SBATCH -J QC
#SBATCH -o /home/twinkle1/projects/Ahyp_v2_2/logs/flower_color_mapping/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=8gb

# Use featurecounts to create the countmatrix from STAR mappings, also perform different quality control measures

# load necessary modules

source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate featurecounts

module load fastqc/0.11.9

GTFIN=polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.gtf



# Quality control
# Initialize the output directory
QCIN=raw_data/BSA_rnaseq/
QCOUT=raw_data/BSA_rnaseq/fastqc

mkdir -p $QCOUT
# run quality control fastqc
fastqc -t 10 -o $QCOUT "$QCIN"*P.fq.gz


# Quality control
# Initialize the output directory
QCIN=raw_data/BSA_rnaseq/
QCOUT=raw_data/BSA_rnaseq/fastqc

mkdir -p $QCOUT
# run quality control fastqc
fastqc -t 10 -o $QCOUT "$QCIN"*P.fq.gz

cp -r $QCOUT data/BSA/RNAseq/STAR_flower_mappings/QC/

module load samtools/1.13


# run quality control qualimap on the generated bam files
# run rnaseq mode for each file, qualimap takes as input a bam sorted by name
# -p = strand specific protocol, -pe = paired-end sequencing data, -s = file is sorted by name

# prepare input gtf file
GTFQM=/scratch/twinkle1/temp.gtf
sed 's/CDS/exon/' $GTFIN > $GTFQM

# define input files
QMIN=data/BSA/RNAseq/STAR_flower_mappings/


# define output directory
QMOUT=data/BSA/RNAseq/STAR_flower_mappings/QC/qualimap

# create main output directory
mkdir -p $QMOUT

# run qualimap
for file in "$QMIN"*out.bam
do
	# basename of each sample
	QMBASE="$(basename -s .sortedByCoord.out.bam $file)"
	# make output directory
	mkdir "$QMOUT"/"$QMBASE"
	# sort by name for qualimap
	samtools sort -n -T /scratch/twinkle1/ -@ 8 $file -O bam > "$QMOUT"/"$QMBASE"/"$QMBASE".name_sorted.bam

	# run qualimap
	qualimap rnaseq -bam "$QMOUT"/"$QMBASE"/"$QMBASE".name_sorted.bam \
		-outdir "$QMOUT"/"$QMBASE" \
		-gtf $GTFQM \
		-p strand-specific-reverse \
		-pe \
		-s \
		--java-mem-size=4G
done

rm $GTFQM

# run multiqc to combine the results from fastqc and qualimap into a single report
MULTIQCOUT=data/BSA/RNAseq/STAR_flower_mappings/multiqc
MULTIQCIN=data/BSA/RNAseq/STAR_flower_mappings/QC/

mkdir -p $MULTIQCOUT

multiqc -o $MULTIQCOUT $MULTIQCIN
