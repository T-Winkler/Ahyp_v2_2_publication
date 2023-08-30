## Merge of computational annotation and Iso-Seq transcripts

Combine the computational annotation with full-length transcript sequencing data using TSEBRA.

### Script order:

- code/merge_annotation/Braker2_subsets_and_merge.Rmd
Prepare BRAKER2 input, use only predicted genes supported by external evidence for the merge

- code/merge_annotation/reannotation_correction.Rmd
Merge BRAKER2 and full-length transcript sequencing data using TSEBRA, deduplicate, rename genes and compare annotation completeness using BUSCO
