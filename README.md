# Banzai Dada2 Pipeline

Bioinformatic pipeline for processing amplicon data using Dada2. Initially based on the Banzai pipeline https://github.com/jimmyodonnell/banzai, this shell script has been adapted at MBARI to use a different set of programs as well as a different taxonomic annotation strategy.

### Input Files:
- Metadata File (.csv)
- Parameter File
- Fastq files of forward and reverse read amplicon sequence data, organized by sample in folders

### Steps:
1. Removes primers using atropos (https://github.com/jdidion/atropos)
2. Runs Dada2 for quality filtering, denoising, chimera removal, read merging and ASV creation (https://benjjneb.github.io/dada2/)
3. Runs Blastn search and uses MEGAN to filter BLAST results, assigning taxonomy (https://github.com/husonlab/megan-ce)
4. Provides further curated products and allows filtering of taxonomic species and genus annotations by bitscore and percent identity

### Ouput Folder:
- Analysis folder named by datetime of run which contains all intermediate and final products as well as copies of the initial input files used to intiate the run.
