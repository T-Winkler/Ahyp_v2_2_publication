#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -o /home/twinkle1/projects/Ahyp_v2_2/logs/bsa/mappingLog-%j.txt
#SBATCH -t 11-00:00:00
#SBATCH -J map_reads
#SBATCH --array=0-9
#SBATCH --nodes=1-1
#SBATCH --ntasks 8
#SBATCH --mem 48g

####SLURM_ARRAY_TASK_ID=0



module use /opt/rrzk/modules/experimental
module load bwamem2/2.0_gnu
module load samtools/1.13


REFERENCE=polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta
#bwa-mem2 index $REFERENCE


PROVIDER=CCG


INPUTPATH=raw_data/BSA_wgs/ #
OUTPUTPATH=data/BSA/wgs/bam_files/ # Change this for different generations

mkdir -p $OUTPUTPATH
mkdir -p ${OUTPUTPATH}/metrics/

#readarray -t FASTQFILESR1 < ${INPUTPATH}/R1read_list.txt
#readarray -t FASTQFILESR2 < ${INPUTPATH}/R2read_list.txt

FASTQFILESR1=($(ls -d $INPUTPATH/*_1.fq.gz))
FASTQFILESR2=($(ls -d $INPUTPATH/*_2.fq.gz))


INFILE_R1="${FASTQFILESR1[$SLURM_ARRAY_TASK_ID]}" #INFILE_R1=${INPUTPATH}/"${FASTQFILESR1[0]}"
INFILE_R2="${FASTQFILESR2[$SLURM_ARRAY_TASK_ID]}" #INFILE_R2=${INPUTPATH}/"${FASTQFILESR2[0]}"

echo $INFILE_R1
echo $INFILE_R2
INDNAME=$(basename $INFILE_R1 _1.fq.gz)

#INDNAME=synDH_$(basename $INFILE |awk -F_ '{print $2}')
#INDNAME=$(basename "$INFILE_R1" .fastq.gz|awk -F_ 'BEGIN{OFS="_";};{print $1;}') #For DH lines with UNIMO in their name ->fastqList2.txt


echo maping reads of $INDNAME
SORTED_NAME=${OUTPUTPATH}/sorted_${INDNAME}.bam
echo $SORTED_NAME


bwa-mem2 mem -t 8 -R '@RG\tID:'${INDNAME}'\tSM:'${INDNAME}'\tCN:'${PROVIDER}'\tPL:illumina' $REFERENCE $INFILE_R1 $INFILE_R2 | samtools sort -O bam -o ${SORTED_NAME}


echo mark duplicates
DEDUP_NAME=${OUTPUTPATH}/${INDNAME}.bam
METRICS_FILE=${OUTPUTPATH}/metrics/${INDNAME}.txt
java -Xmx40g -jar /projects/mstette2/tools/picard.jar MarkDuplicates INPUT=${SORTED_NAME} OUTPUT=${DEDUP_NAME} METRICS_FILE=${METRICS_FILE}
samtools index $DEDUP_NAME

echo calculate samtools flagstat
samtools flagstat $DEDUP_NAME > ${OUTPUTPATH}/metrics/${INDNAME}.flagstat

echo removing sorted bam
#rm $SORTED_NAME


mkdir -p ${OUTPUTPATH}/gvcf
GVCFFILE=${OUTPUTPATH}/gvcf/${INDNAME}.g.vcf

$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx16G" CreateSequenceDictionary -R $REFERENCE

$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" HaplotypeCaller \
-R $REFERENCE \
-I $DEDUP_NAME \
-ERC GVCF \
--output ${GVCFFILE}
