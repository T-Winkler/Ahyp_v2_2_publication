## Annotation analysis

Identify betalain and flavonoid pathway genes, as well as MYB transcription factor genes in the new genome annotation v2.2.

### Script order:

- code/annotation_analysis/circos_plotting.Rmd
Generate circos plot based on gene and repetitive element annotation

- code/annotation_analysis/betalain_and_flavonoid_identification.sh
Identify candidate betalain genes by BLAST, candidate flavonoid genes by KIPEs

- code/annotation_analysis/betalain_phylogenetic_analysis.Rmd
Phylogenetic analysis of betalain pathway genes

- code/annotation_analysis/R2R3_analysis_reannotation.Rmd
Identify MYB transcription factor genes by HMMscan, analyse putative function by phylogenetic assessment
