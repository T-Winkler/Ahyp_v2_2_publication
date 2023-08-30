#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 01:00:00
#SBATCH -J minimap
#SBATCH -o /home/twinkle1/projects/Ahyp_v2_2/logs/other_inferences/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=20gb


# load necessary modules
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate isoseq

# create output directory
MAPPINGOUT=data/mapping_bias_inference/plastid_to_genome/
mkdir -p "$MAPPINGOUT"

### MAPPING
# map reads onto the reference genome
REFERENCE=/projects/ag-stetter/reference_genomes/Ahypochondriacus/V2_2/Ahypochondriacus_2.2_polished.softmasked.fasta

# Align sequences to reference genome
# sequences obtained from (Beta vulgaris, mitchondrium): https://www.ncbi.nlm.nih.gov/nuccore/BA000009.3
# (Amaranthus hypochondriacus, chloroplast): https://www.ncbi.nlm.nih.gov/nuccore/KX279888.1

INPUTMITO="$MAPPINGOUT"Bv_mitochondrium.fasta
INPUTCHLORO="$MAPPINGOUT"Ah_chloroplast.fasta
OUTPUTMITO="$MAPPINGOUT"mitochondrium_to_genome.paf
OUTPUTCHLORO="$MAPPINGOUT"chloroplast_to_genome.paf


# asm5 for intra species/genus chloroplast assembly, asm10 for cross-species mitchondrial alignment
minimap2 -t 8 -cx asm10 $REFERENCE $INPUTMITO > $OUTPUTMITO
minimap2 -t 8 -cx asm5 $REFERENCE $INPUTCHLORO > $OUTPUTCHLORO
