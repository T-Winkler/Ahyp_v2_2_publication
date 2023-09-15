#!/bin/bash -l
#SBATCH -D /projects/ag-stetter/twinkle/projects/Ahyp_v2_2_publication/
#SBATCH -o /projects/ag-stetter/markus/bsa_sterility_color/logs/callingLog-%j.txt
#SBATCH -t 14:00:00
#SBATCH -J map_reads
#SBATCH --nodes=1-1
#SBATCH --ntasks 5
#SBATCH --mem 48g

#module load bwa
#module load java/1.8
module load samtools/1.13


REFERENCE=polished_genome_annotation/assembly/Ahypochondriacus_2.2_polished.softmasked.fasta

PROVIDER=CCG
OUTDIR=data/BSA/wgs/vcf
mkdir -p $OUTDIR

ALLSAMP=$(for i in data/BSA/wgs/bam_files/gvcf/AM_00*.vcf; do echo -V $i;done)

$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" \
	CombineGVCFs \
   -R $REFERENCE \
   $ALLSAMP \
   -O $OUTDIR/cohort.g.vcf.gz


$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" \
GenotypeGVCFs \
-R $REFERENCE \
-V $OUTDIR/cohort.g.vcf.gz \
-O $OUTDIR/raw_snps_all.g.vcf \
--sample-ploidy 50  # this is for pool data pools are approx 25 ind


$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" \
VariantFiltration \
-R $REFERENCE \
-V $OUTDIR/raw_snps_all.g.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_snp_filter" \
--output $OUTDIR/raw_variants_gatk.vcf

$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" \
SelectVariants \
-R $REFERENCE \
-V $OUTDIR/raw_variants_gatk.vcf  \
--select-type-to-include SNP \
--output $OUTDIR/filtered_snps_gatk.vcf

vcftools --vcf $OUTDIR/filtered_snps_gatk.vcf \
--remove-filtered-all --min-alleles 2 --max-alleles 2 --max-missing 0.95 --recode \
--out $OUTDIR/gatk_filter_maxmissing05_biallelic

mv $OUTDIR/gatk_filter_maxmissing05_biallelic.recode.vcf $OUTDIR/gatk_filter_maxmissing05_biallelic.vcf

$MYUTIL/tools/gatk-4.1.7.0/gatk --java-options "-Xmx40G" \
VariantsToTable \
-R $REFERENCE \
-V $OUTDIR/gatk_filter_maxmissing05_biallelic.vcf \
-F CHROM -F POS -F REF -F ALT \
-GF AD -GF DP -GF GQ -GF PL \
--output $OUTDIR/bulk_snps05.table
