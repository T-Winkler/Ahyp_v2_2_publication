#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 1:00:00
#SBATCH -J STAR
#SBATCH -o logs/STAR/braker/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --job-name="index_STAR"

module load star/2.7.8a

# Index the reference genome
# only run once per reference genome

# 8 threads, genome generation mode
# sjdbOverhang and genomeSAindexNbases settings specific for the amaranth reference assembly v2.1
# more specific settings: use the polished, softmasked reference assembly
# as SJDB file, use the newly generated braker2 protein gtf file


mkdir -p data/braker2/STAR_index/

STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir data/braker2/STAR_index/ \
	--sjdbOverhang 89 \
	--genomeSAindexNbases 13 \
	--genomeFastaFiles polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta \
	--sjdbGTFfile data/braker2/prot_run/braker.gtf
