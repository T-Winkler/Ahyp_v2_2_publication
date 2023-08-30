#!/bin/bash -l
#SBATCH -D /scratch/twinkle1/nextpolish/
#SBATCH -t 40:00:00
#SBATCH -J nextpolish
#SBATCH -o /home/twinkle1/master_thesis/logs/nextpolish/mappingLog-%j.txt
#SBATCH --error /home/twinkle1/master_thesis/logs/nextpolish/errorLog-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=21
#SBATCH --mem=100gb


# Run nextpolish to error correct a reference assembly using short or long reads. Adjust parameters under "Set input and parameters"


# setup
module load bwamem2/2.2.1
module load samtools

mkdir -p /scratch/twinkle1/nextpolish/

# set output directory for saving the polished genome to:
OUTDIR=data/NextPolish/output/



### Prepare input files:
# remove all reads containing ambiguous bases from the input using bbduk
/home/twinkle1/tools/bbmap/bbduk.sh maxns=0 \
	in=/projects/ag-stetter/twinkle/lightfoot_WGS_short_reads/SRR2106212/SRR2106212_1.fastq.gz \
	in2=/projects/ag-stetter/twinkle/lightfoot_WGS_short_reads/SRR2106212/SRR2106212_2.fastq.gz \
	out=/scratch/twinkle1/SRR2106212_1.cleaned.fq \
	out2=/scratch/twinkle1/SRR2106212_2.cleaned.fq \
	-Xmx16g

# make sure all reads are still paired and the files are in the correct order by running repair from the bbmap suite:
/home/twinkle1/tools/bbmap/repair.sh \
	in=/scratch/twinkle1/SRR2106212_1.cleaned.fq \
	in2=/scratch/twinkle1/SRR2106212_2.cleaned.fq \
	out=/scratch/twinkle1/SRR2106212_1.cleaned.repair.fq \
	out2=/scratch/twinkle1/SRR2106212_2.cleaned.repair.fq \
	-Xmx16g



### Adopted script to use NextPolish manual from: https://nextpolish.readthedocs.io/en/latest/TUTORIAL.html

#Set input and parameters
round=2
threads=20
read1=/scratch/twinkle1/SRR2106212_1.cleaned.repair.fq
read2=/scratch/twinkle1/SRR2106212_2.cleaned.repair.fq
input=/home/twinkle1/master_thesis/data/NextPolish/input/Ahypochondriacus_split.fasta


for ((i=1; i<=${round};i++)); do
#step 1:
   #index the genome file and do alignment
   bwa-mem2 index ${input};
   bwa-mem2 mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 19 -F 0x4 -b -|samtools fixmate -m --threads 19  - -|samtools sort -m 2g --threads 20 -|samtools markdup --threads 19 -r - sgs.sort.bam
   #index bam and genome files
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${input};
   #polish genome file
   python /home/twinkle1/tools/NextPolish/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam -debug > genome.polishtemp.fa;
   input=genome.polishtemp.fa;
#step2:
   #index genome file and do alignment
   bwa-mem2 index ${input};
   bwa-mem2 mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 19 -F 0x4 -b -|samtools fixmate -m --threads 19  - -|samtools sort -m 2g --threads 20 -|samtools markdup --threads 19 -r - sgs.sort.bam
   #index bam and genome files
   samtools index -@ ${threads} sgs.sort.bam;
   samtools faidx ${input};
   #polish genome file
   python /home/twinkle1/tools/NextPolish/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam -debug > genome.nextpolish.fa;
   input=genome.nextpolish.fa;
done;
#Finally polished genome file: genome.nextpolish.fa

cp /scratch/twinkle1/nextpolish/* $OUTDIR
