## Bulk segregant analysis

Code for bulk segregant analysis, haplotype analysis and quantification of flower gene expression.


### Script order

RNAseq:

- code/BSA/RNAseq/adapter_trimming.sh
trimming of adapter sequences using Trimmomatic

- code/BSA/RNAseq/index_STAR.sh
index the genome with STAR using the genome annotation v2.2

- code/BSA/RNAseq/run_STAR.sh
map reads to the genome using STAR

- code/BSA/RNAseq/QC.sh
quality control of adapter trimming and read mapping

- code/BSA/RNAseq/run_kallisto.sh
quantification of gene expression using kallisto

WGS:

- code/BSA/WGS/map_reads.sh
map WGS reads to the genome

- code/BSA/WGS/combined_filter.sh
call and filter variants

Analysis:

- code/BSA/read_count_analysis.Rmd
Analysis of gene expression data in pooled flower tissue

- code/BSA/bsa_and_plotting.R
Conduct bulk segregant analysis

- code/BSA/snpEff_analysis.Rmd
Analyse BSA variants using snpEff

- code/BSA/read_count_analysis_from_bam.Rmd
Investigate support for AhCYP76AD2 variants by RNA-seq data
