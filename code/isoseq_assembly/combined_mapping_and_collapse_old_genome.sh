#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -t 02:00:00
#SBATCH -J collapse
#SBATCH -o logs/mapping_and_collapse/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=42gb


# load necessary modules
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate isoseq

# create output directory
MAPPINGOUT=data/isoseq/mapping_and_collapse_old_genome/
mkdir -p "$MAPPINGOUT"


### MAPPING
# map reads onto the old reference genome
REFERENCE=/home/tom/Documents/reference_genomes/Ahypochondriacus/assembly/Ahypochondriacus_459_v2.0.softmasked.nospace.underscore.fa

# Align sequences to reference genome
INPUT=raw_data/isoseq_raw_reads/processed/combined.merged_clustered.hq.fasta
OUTPUTMM2="$MAPPINGOUT"combined_aln.sam

minimap2 -t 8 -ax splice:hq -uf --secondary=no -a $REFERENCE $INPUT -o $OUTPUTMM2



# Before collapsing isoforms, sequences have to be sorted
OUTPUTSORT="$MAPPINGOUT"combined_aln_sorted.sam

# remove unmapped sequences from file and sort
sort -k 3,3 -k 4,4n $OUTPUTMM2 > $OUTPUTSORT



### COLLAPSE

# activate conda environment
conda activate /projects/ag-stetter/twinkle/sqanti_env/sqanti


# Collapsing isoforms

INPUTSORTED="$MAPPINGOUT"combined_aln_sorted.sam
OUTPUTCOLLAPSE="$MAPPINGOUT"combined

collapse_isoforms_by_sam.py --input $INPUT -s $INPUTSORTED -o $OUTPUTCOLLAPSE -c 0.95 -i 0.9 --max_3_diff 1000



# Cupcake support scripts after collapse
# First obtain associated count information

INPUTABUNDANCE="$MAPPINGOUT"combined.collapsed
CLUSTERREPORT=raw_data/isoseq_raw_reads/processed/combined.merged_clustered.cluster_report.csv

get_abundance_post_collapse.py $INPUTABUNDANCE $CLUSTERREPORT

# add minimum read count of 2

filter_by_count.py --min_count 2 --dun_use_group_count $INPUTABUNDANCE


# filter away 5' degraded isoforms
# use filter by count results

OUTPUTFILTERED="$MAPPINGOUT"combined.collapsed.min_fl_2

filter_away_subset.py $OUTPUTFILTERED
