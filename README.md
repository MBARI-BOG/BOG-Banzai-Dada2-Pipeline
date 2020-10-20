# Banzai Dada2 Pipeline

Bioinformatic pipeline for processing amplicon data using Dada2. Initially based on the Banzai pipeline https://github.com/jimmyodonnell/banzai, this shell script has been adapted at MBARI to use a different set of programs as well as a different taxonomic annotation strategy.

### Input Files:
- Metadata File (.csv)
- Parameter File
- Fastq files of forward and reverse read amplicon sequence data, organized by sample in folders

### Overview:
1. Primer sequences are removed from fastq files through atropos (https://github.com/jdidion/atropos)
2. Fastq files are then fed into the Dada2 program for quality filtering, denoising, chimera removal, read merging and ASV creation (https://benjjneb.github.io/dada2/)
3. Taxonomy is assigned through blastn searches to NCBI GenBank’s non-redundant nucleotide database (nt) and hits are filtered by MEGAN6 and custom scripts, as described below (https://github.com/husonlab/megan-ce)

### Ouput Folder:
- Analysis folder named by datetime of run which contains all intermediate and final products as well as copies of the initial input files used to intiate the run.

### Detailed Description

The Banzai Dada2 Pipeline is a custom shell script adapted from the banzai pipeline and run on a per-plate basis superscript^1^. Within the script, primer sequences are first removed from fastq files through the program Atropos superscript^2^. Fastq files are then fed into the Dada2 program superscript^3^. Dada2 models error on a per-Illumina run basis, controlling for read quality and picking amplicon sequence variant (ASV) sequences that represent biological variability rather than sequencing error superscript^3^. Within Dada2 reads are trimmed to remove low-quality regions and filtered by quality score, sequencing errors are modeled and removed, and reads are then merged and chimeric sequences removed. Taxonomy is assigned to the resulting ASV sequences through blastn searches to NCBI GenBank’s non-redundant nucleotide database (nt) superscript^4,5^. Blast results are filtered using MEGAN6’s lowest common ancestor (LCA) algorithm superscript^6^. Only hits with ≥80% sequence identity and whose bitscores were within the top 2% of the highest bitscore value for each OTU are considered by MEGAN6. The MEGAN6 parameter LCA percent is set to 0.8, allowing for up to 20% of top hits to be off target and still have the majority taxonomy assigned. This parameter value was chosen to allow for minor numbers of incorrectly annotated GenBank entries – effectively allowing for OTUs which had many high-quality hits to a taxa to still be assigned to that taxa even if there existed a high-bitscore hit to another GenBank sequence annotated to an unrelated taxa. We decided this was more advantageous than the disadvantage caused by ignoring small numbers of true closely related sequences. Furthermore, post-MEGAN6 filtering is performed to ensure only contigs with a hit of ≥97% sequence identity are annotated to the species level. Only contigs with a hit of ≥95% sequence identity are annotated to the genus level. Annotations are elevated to the next highest taxonomic level for contigs that fail these conditions.

1. 	O’Donnell JL, Kelly RP, Lowell NC, Port JA. Indexed PCR primers induce template- Specific bias in Large-Scale DNA sequencing studies. PLoS One. 2016;11(3):1–11.
2. 	Didion JP, Martin M, Collins FS. Atropos: Specific, sensitive, and speedy trimming of sequencing reads. PeerJ. 2017 Aug 30;2017(8):e3720.
3. 	Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. 2016 Jun 29;13(7):581–3.
4. 	Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: architecture and applications. BMC Bioinformatics. 2009 Dec 15;10(1):421.
5. 	Agarwala R, Barrett T, Beck J, Benson DA, Bollin C, Bolton E, et al. Database resources of the National Center for Biotechnology Information. Nucleic Acids Res. 2018 Jan 1;46(D1):D8–13.
6. 	Huson DH, Beier S, Flade I, Górska A, El-Hadidi M, Mitra S, et al. MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data. PLoS Comput Biol. 2016;12(6):1–12.
