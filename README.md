# Computational_Biology
 This is normal text.

## Snakemake Workflows
These are designed to run on the NIH HPC (Biowulf set-up using SLURM). With a few modifications to the shell scripts and to the module commands (e.g. loading sratoolkit) these can be run on another system.

- *star_rsem_L75_SE* -- Takes as input SRA acc ids for single-end 75-base RNA-seq files (Illumina) and uses the STAR-RSEM pipeline to output RNA expression quantification at gene- and isoform-levels for each sample. Quality control performed with Trimmomatic.
- *MHC_query_predict* -- Takes as input a headerless list (single column) of peptides and uses pvacbind to output HLA allotype-specific binding predictions using NetMHC and MHCflurry for each peptide length in the provided range.
- *standard_MQ_grid* -- Takes as input FASTA protein databases, mass spectrometry RAW files, and mqpar.xml configuration files for MaxQuant. Runs MaxQuant using each combination of protein database, peptide size range, and MaxQuant FDR style, and outputs MaxQuant results.

## Analyses in R
These can be run locally in an IDE like *RStudio*
- *GSEA_tutorial* -- Gene Set Enrichment Analysis tutorial using *R* package *fgsea*
- *Survival_analysis* -- *R* script running parallel multivariate multilevel Cox proportional hazards model analyses using TCGA SKCM gene- and isoform-level RNA-seq expression data and TCGA clinical data. (Script is part of a larger workflow.)

## Utilities in C++
These can be run in the command line once they are compiled on the system in use
- *build_kmers_simple.cpp* -- Build a list of unique k-mers based on a converted FASTA sequence file
- *build_IL_tbl.cpp* -- Convert a list of sequences to matched pairs of the original sequence and the sequence with every 'I' replaced with 'L'

- recursive letter swaps
- in silico translate NT to AA
- search for string in protein file

## Utilities in AWK
These can be run in the command line
- *FastaToTbl.awk* -- Converts a FASTA file to a .tbl file (converts every two lines in the input FASTA file to a single two-column tab-separated line in the output .tbl file); copied from a StackExchange post: https://bioinformatics.stackexchange.com/questions/2649/how-to-convert-fasta-file-to-tab-delimited-file*
- *TblToFasta.awk* -- Converts a .tbl file to a FASTA file (converts each two-column tab-separated line in the input .tbl file to two lines where the first starts with ">"); copied from a StackExchange post: https://bioinformatics.stackexchange.com/questions/2649/how-to-convert-fasta-file-to-tab-delimited-file