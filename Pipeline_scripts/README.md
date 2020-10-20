# Banzai Dada2 Pipeline Scripts

## Overview:

- **banzai_dada2_2019v1.4.sh** : overall bash script that links together bioinformatic programs and other scripts
- **Dada2_2019v1.0.R** : Runs Dada2 program, called within **banzai_dada2_2019v1.4.sh**
- **blastn-wrapper** : part of the virtual machine BLAST strategy - script resides on "worker" VMs
- **Make_otu_taxa_table.py** : runs within **banzai_dada2_2019v1.4.sh**, takes output of Dada2 and MEGAN and formats into csv files to create the ASV and Taxa tables.
- **XML_Filter_dada2.py** : Creates optional "Filtered" ASV and Taxa tables where ASVs with hits to PhiX are removed and imposes limits on bitscore and Percent Identity for ASV assignment at the species or genus level (set within the input parameter file)
- **Post_BLAST_processing_dada2_v2.0.sh** : bash script linking together scripts and programs for processing post-blastn products. This script is copied to the Analysis folder alongside parameters specific to each run after the blastn search concludes. The user then runs this script to finish the pipeline.
- **decon_std_outputs_v1.0.R** : An optional script which takes the outputs of the bioinformatic pipeline and runs through decontamination steps, removes singletons, and produces some standard plots. Is called within **Post_BLAST_processing_dada2_v2.0.sh**.


## Virtual Machine BLAST strategy

In order to speed up BLAST taxonomic annotations, ASVs generated through Dada2 are split into 6 input files that are run through GNU parallel to 6 "worker" virtual machines which run the blastn search in parallel. This processes is controlled by the scripts banzai_dada2_2019v1.4.sh and blastn-wrapper.

Since the blastn search is time intensive, this pipeline is split into Pre- and Post- BLAST bash analysis scripts, enabling re-running of Post-BLAST analyses without having to repeat the initial Dada2 and BLAST processes.
