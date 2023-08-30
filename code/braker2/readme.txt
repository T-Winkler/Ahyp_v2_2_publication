## Computational annotation

Create computational genome annotation using BRAKER2.
 
### Script order:

- code/braker2/braker2_prot.sh
initial run of BRAKER2 using only the protein database as evidence

- code/braker2/index_STAR.sh
index the polished reference genome for use with STAR

- code/braker2/run_STAR.sh
map all short reads from Clouse et al. to the polished reference genome for use with BRAKER2

- code/braker2/braker2_prot_rna.sh
run BRAKER2 with the mapped RNA-seq reads as well as the protein database as input
