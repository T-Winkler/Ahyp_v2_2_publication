## Iso-Seq assembly

Assembly full-length transcript sequencing data.

### Script order:

- code/isoseq_assembly/isoseq3_pipeline.sh
assemble FLNC reads from CCS files

- code/isoseq_assembly/combined_isoseq3_pipeline.sh
combine FLNC reads and cluster identical reads

- code/isoseq_assembly/combined_mapping_and_collapse.sh
collapse clustered reads into unique full-length transcripts using the polished reference genome

- code/isoseq_assembly/combined_mapping_and_collapse_old_genome.sh
collapse clustered reads into unique full-length transcripts using the unpolished reference genome

- code/isoseq_assembly/run_sqanti.sh
run SQANTI in order to correct possible sequencing errors in the full-length transcripts using both reference genomes

- code/isoseq_assembly/comparison_isoseq_polishing_effectiveness.Rmd
compare effect of genome polishing on error correction of full-length transcript sequences
