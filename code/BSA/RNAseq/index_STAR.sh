#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 1:00:00
#SBATCH -J STAR
#SBATCH -o /home/twinkle1/projects/Ahyp_v2_2/logs/flower_color_mapping/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --job-name="index_STAR"

module load star/2.7.8a

# Index the reference genome
# only run once per reference genome

# 8 threads, genome generation mode
# genomeSAindexNbases settings specific for the amaranth reference assembly
# more specific settings: use the polished, softmasked reference assembly
# sjdbOverhang dependend on input read length
# as SJDB file, use the newly generated braker2 protein gtf file
REFGENOME=polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta
SJDBFILE=polished_genome_annotation/annotation/Ahypochondriacus_2.2_polished_corrected.gtf

mkdir -p data/BSA/RNAseq/STAR_flower_index

STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir /scratch/twinkle1/STAR_flower_index \
	--sjdbOverhang 100 \
	--genomeSAindexNbases 13 \
	--genomeFastaFiles "$REFGENOME" \
	--sjdbGTFfeatureExon CDS \
	--sjdbGTFfile "$SJDBFILE"
