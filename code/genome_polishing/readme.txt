## Genome polishing

Polish the previously published reference genome of A. hypochondriacus.

### Script order:

- code/genome_polishing/remove_ambiguous_bases.sh
prepare the reference genome v2.1 for genome polishing by removing all ambiguous bases

- code/genome_polishing/unpackSRA.sh 
unpack the WGS short read SRA file from the Lightfoot genome assembly

- code/genome_polishing/run_nextpolish.sh
repare input files and run NextPolish to polish the reference assembly.

- code/genome_polishing/process_nextpolish_output.sh
return all ambigious bases into the polished reference genome and restore the same chromosome order


The helper_script.R is called by other scripts and does not have to be manually run.
