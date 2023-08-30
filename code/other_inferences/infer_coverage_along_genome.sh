#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2
#SBATCH -t 30:00:00
#SBATCH -J genome_coverage
#SBATCH -o logs/other_inferences/mappingLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=48gb


# Run bwa to map WGS to the reference genome, calculate read depth at each position using samtools

# setup
module load bwamem2/2.2.1
module load samtools/1.13


#Set input and parameters
read1=raw_data/lightfoot_WGS_short_reads/SRR2106212/SRR2106212_1.fastq.gz
read2=raw_data/lightfoot_WGS_short_reads/SRR2106212/SRR2106212_2.fastq.gz
input=polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta
outdir=data/mapping_bias_inference/genome_coverage/
outbam=SRR2106212.sorted.bam
outdedup=SRR2106212.sorted.dedup.bam
outcov=SRR2106212.coverage

mkdir -p outdir

#index the genome file and map
bwa-mem2 index ${input}
# map using 10 threads, output in bam format with header, fixmates, sort and mark and remove duplicates
bwa-mem2 mem -t 10 ${input} ${read1} ${read2}  | samtools sort -O bam -o "$outdir""$outbam"

echo "mapping complete"

# load correct java version
module load openjdk/1.8.0_60

echo mark duplicates

java -Xmx40g -jar /home/twinkle1/tools/picard/picard.jar MarkDuplicates --TMP_DIR /scratch/twinkle1 --INPUT "$outdir""$outbam" --OUTPUT "$outdir""$outdedup" --METRICS_FILE "$outdir"SRR2106212.sorted.dedup.metrics


samtools index -@ 9 "$outdir""$outdedup"

echo "marking duplicates complete"

#index bam and genome files
#samtools index -@ 9 $outbam

# get flagstat report
samtools flagstat -@ 9 "$outdir""$outbam" > "$outdir""$outbam".flagstat.txt

echo "indexing complete"

# calculate depth
samtools depth "$outdir""$outdedup" > "$outdir"SRR2106212.sorted.dedup.depth

# subset scaffold 10
grep "Scaffold_10" "$outdir"SRR2106212.sorted.dedup.depth > "$outdir"SRR2106212.sorted.dedup.10.depth

echo "finished"
