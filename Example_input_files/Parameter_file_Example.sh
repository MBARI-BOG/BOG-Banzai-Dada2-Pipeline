#!/usr/bin/env bash


################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?
PARENT_DIR="/MBON/raw_sequence_data/18S/S_20170403_18S_PE_M2W/"

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="/home/mbonteam/MBARI/kpitz/scripts/S_metadata.csv"

#Value can be either "COI","12S", or "18S"
MARKER="12S"
################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
# Create this directory before running Banzai
ANALYSIS_DIRECTORY="/home/mbonteam/MBARI/kpitz/processed/S"

# You can optionally specify a folder into which the script copies a PDF containing some results.
# The pdf is created by default in the analysis folder specified above, but
# if you set this to your DropBox or Google Drive Folder, you can check it out from anywhere.
OUTPUT_PDF_DIR="/home/mbonteam/MBARI/kpitz/processed/S"
SANCTUARY="MB"
LOCUS="18S"
BIOM_FILE_NAME="S_18S"

################################################################################
# METADATA DETAILS
################################################################################
#Expected length of reads after trimming by cutadapt (not full product, but read)
#18S: 120  #150 seq, 18 and 16nt primers - BUT SEQUENCE FACILITY SHOULD REMOVE PRIMERS
#COI: 224  #250 seq, 26nt primers
#12S: 171   #250 seq, ~171 nt Euk final length
#READ_LEN="171"   #In the future can adjust cutadapt min/max based on read length

# Your metadata must have a column corresponding to the subfolders containing the raw reads.
# In order to make this flexible across both multiple and single library preps, you must include this even if you only sequenced one library (sorry!).
READ_LIB_FROM_SEQUENCING_METADATA="YES"
LIBRARY_COLUMN_NAME="library"
LIBRARY_TAG_COMBO_COLUMN_NAME="library_tag_combo"

################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex tags, Banzai will grab these from the file SEQUENCING_METADATA.
# You must indicate the column names of the forward and reverse primers
REMOVE_PRIMER="YES"
PRIMER_1_COLUMN_NAME="primer_sequence_F"
PRIMER_2_COLUMN_NAME="primer_sequence_R"

# What proportion of mismatches are you willing to accept when looking for primers?
# Recommended: "0.10"
#PRIMER_MISMATCH_PROPORTION="0.10"

ColumnName_SampleName="sample_name"
ColumnName_SampleType="sample_type"

################################################################################
# PRIMER REMOVAL
################################################################################
#Dada2

#Based on read quality scores, set primer truncation lengths for dada2
#typical values:
#COI F- 200, R - 200
#12S F - 160,R - 160  #euk len without primers ~ 173
#18S F - 115, R - 115

truncF="115"
truncR="115"

################################################################################
# TAXONOMIC ANNOTATION
################################################################################
## BLAST ##
# For more information on these parameters, type into a terminal window: blastn -help
# Specify the path to the BLAST database.
# Note this should be a path to any one of three files WITHOUT their extension *.nhr, *.nin, or *.nsq
#BLAST_DB='/MBON/blastdb/nr/nr'
#BLAST_DB='/MBON/blastdb/nr/nr'
BLAST_DB='/MBON/blastdb/nt/nt'
#BLAST_DB='/MBON/blastdb/Tarav9/Tarav9'
#BLAST_DB='/MBON/blastdb/greengenes/gg_13_5'
#BLAST_DB='/MBON/blastdb/greengenes/gg_13_5_with_header'
#BLAST_DB='/MBON/blastdb/Silva/SILVA_128_SSURef'
# BLAST PARAMETERS
PERCENT_IDENTITY="80"
#18S WORD_SIZE = 20
#COI WORD_SIZE = 11
WORD_SIZE="11"
EVALUE="1e-5"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="100"
#culling_limit="20"   #REMOVED from search by BLAST recommendations KP
#best_hit_score_edge="0.05"
#best_hit_overhang="0.25"


# ::::::::::::::::::: Filter XML File for PhiX and apply genus, species limits ::::::::::::::::::::::::
#Filter out OTUs with hits to PhiX - remove from complete OTU table (will be unannotated)
#Important for post blast processing script
#Also need to filter annotations at the genus and species level to different bitscore, %ID limits

#BLAST parameters for species assignment:
#18S bitscore ~ 200 ; COI bitscore ~ 400 (COI has a longer alignment length)
bitscore_sp=200
per_ID_sp=97
#BLAST parameters for genus assignment:
#18S bitscore ~ 150 ; COI bitscore ~ 400 (COI has a longer alignment length)
bitscore_gn=150
per_ID_gn=95

################################################################################
## MEGAN ##   #This section may be important for the Post blast processing script?
# For more information, see the manual provided with the software
# Specify the path to the MEGAN executable file you want to use.
# Note that in recent versions an executable was not provided; in that case, you need to reference like so: '/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub'
#megan_exec='/usr/local/bin/MEGAN'

# What is the lowest taxonomic rank at which MEGAN should group OTUs?
# COLLAPSE_RANK1="Family"
 MINIMUM_SUPPORT="1"
# MINIMUM_COMPLEXITY="0"
 TOP_PERCENT="2"
# MINIMUM_SUPPORT_PERCENT="0"
 MINIMUM_SCORE="100"
 LCA_PERCENT="80"
 MAX_EXPECTED="1e-25"


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to compress extraneous intermediate files once the analysis is finished? YES/NO
#PERFORM_CLEANUP="YES"

# If you want to receive a text message when the pipeline finishes, input your number here:
NOTIFY_EMAIL="YES"
EMAIL_ADDRESS="youremail@site.com"
