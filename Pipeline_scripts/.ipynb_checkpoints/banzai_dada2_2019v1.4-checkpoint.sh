#!/usr/bin/env bash
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# banzai_2019v1.2.sh
#
# This script is customized to run the MBARI,STANFORD,USF sequences
# updated Sep 2017
# R. Michisaki MBARI
# CUTADAPT. Major reconstruction to run cutadapt for MBARI samples. Add param INDEXED. Put primer removal, cutadapt
# in library level loop.
# Aug 2019 rpm
# Dada2. Major changes to run Dada2 instead of swarm. KP Nov 19 2019
# Cutadapt replaced with Atropos KP 011620
# Atropos parameters modified to trim primer repeats 022520 KP
# Changed permissions structure so all files generated have permissions open to everyone
# Now calling raw data from /atlasMBON/ exclusively
# Updated September 2020 by Markus Min to now include decontamination and standard outputs as part of the pipeline.
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Pipeline for analysis of MULTIPLEXED Illumina data, a la Jimmy
# Modified for the USF 16S sequences.  Modifications were made by Reiko Michisaki MBARI, Dec 2015, Jan 2016
# Modified to use Dada2 instead of swarm by KP Nov 19 2019

echo
echo
echo -e '\t' "\x20\xf0\x9f\x8f\x84" " "  "\xc2\xa1" BANZAI plus Dada2 !
echo
echo

# Change permissions of all files generated: set to 666 if text or 777 if executable
umask 000
#directory 775, file 664
#umask 002
################################################################################
# CHECK FOR RAW DATA
################################################################################

# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)

# Find the directory this script lives in, so it can find its friends.
#SCRIPT_DIR="$(dirname "$0")"
SCRIPT_DIR="/home/mbonteam/dev"
MY_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"    #where you're running this script?

# Read in the parameter file (was source "$SCRIPT_DIR/banzai_params.sh"; now argument 1)
param_file="${1}"
source "${param_file}"

# check if param file exists:
if [[ -s "${param_file}" ]] ; then
	echo "Reading analysis parameters from:"
	echo "${param_file}"
else
	echo
	echo 'ERROR! Could not find analysis parameter file. You specified the file path:'
	echo
	echo "${param_file}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi


# check if sequencing metadata exists
if [[ -s "${SEQUENCING_METADATA}" ]] ; then
	echo "Reading sequencing metadata from:"
	echo "${SEQUENCING_METADATA}"
else
	echo
	echo 'ERROR! Could not find sequencing metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi
# make standard name prefix for OTU, XML, and biom files
STANDARD_PREFIX="${SANCTUARY}"_"${START_TIME}"_"${LOCUS}"
DB_NAME=$(basename $BLAST_DB)
#BIOM_FILE_NAME="${STANDARD_PREFIX}"_"${DB_NAME}"

# make an analysis directory with starting time timestamp
ANALYSIS_DIR="${ANALYSIS_DIRECTORY}"/Analysis_"${START_TIME}"
mkdir "${ANALYSIS_DIR}"

# Write a log file of output from this script (everything that prints to terminal)
LOGFILE="${ANALYSIS_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1

################################################################################
# Setup Analysis Directory
################################################################################

# Specify compression utility
if hash pigz 2>/dev/null; then
	ZIPPER="pigz"
	echo "pigz installation found"
else
	ZIPPER="gzip"
	echo "pigz installation not found; using gzip"
fi


echo "Analysis started at ""${START_TIME}" " and is located in ""${ANALYSIS_DIR}"

echo "Running script from: "
echo "$BASH_SOURCE"
echo "Parameters read from: "
echo "${param_file}"
echo "metadata read from: "
echo "${SEQUENCING_METADATA}"

# Copy these files into that directory as a verifiable log you can refer back to.
cp $BASH_SOURCE "${ANALYSIS_DIR}"/"${STANDARD_PREFIX}"_analysis_script.txt
cp "${param_file}" "${ANALYSIS_DIR}"/"${STANDARD_PREFIX}"_analysis_parameters.txt
cp "${SEQUENCING_METADATA}"  "${ANALYSIS_DIR}"/"${STANDARD_PREFIX}"_analysis_metadata.csv

################################################################################
# CHECK METADATA FILE FORMAT AND FOR DUPLICATE SAMPLE  NAMES
################################################################################

echo 'Checking for duplicate sample names.'
SID_COL=$(awk -F',' -v SID_COL_NAME=$ColumnName_SampleName '{for (i=1;i<=NF;i++) if($i == SID_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
SID=$(awk -F',' -v SIDCOL=$SID_COL 'NR>1 {print $SIDCOL}' $SEQUENCING_METADATA | uniq)
duplicates=($(printf "${SID}" | sort | uniq -d))

if (( ${#duplicates[@]} )); then
	echo " "
  echo 'FATAL ERROR:' "${ColumnName_SampleName}" 'duplicates found.' "${ColumnName_SampleName}" 'must be unique or the biom conversion will fail.'
	echo 'Duplicate(s):' "${duplicates[@]}"
	echo 'FATAL ERROR: Exiting program'
	echo " "
	exit 1
else
	echo " "
fi

################################################################################
# Read in primers and create reverse complements. R2_COLUMN_NAME="R2"
################################################################################
PRIMER1_COLNUM=$(awk -F',' -v PRIMER1_COL=$PRIMER_1_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == PRIMER1_COL) print i; exit}' $SEQUENCING_METADATA)
PRIMER2_COLNUM=$(awk -F',' -v PRIMER2_COL=$PRIMER_2_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == PRIMER2_COL) print i; exit}' $SEQUENCING_METADATA)
PRIMER1=$(awk -F',' -v PRIMER1_COL=$PRIMER1_COLNUM 'NR==2 {print $PRIMER1_COL}' $SEQUENCING_METADATA)
PRIMER2=$(awk -F',' -v PRIMER2_COL=$PRIMER2_COLNUM 'NR==2 {print $PRIMER2_COL}' $SEQUENCING_METADATA)


echo "Primers read from sequencing metadata:" "${PRIME1}" "${PRIME2}"

if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
  echo 'Primers read from metadata columns' "${PRIMER1_COLNUM}" 'and' "${PRIMER2_COLNUM}"
  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
else
  echo 'ERROR:' 'At least one primer is not valid'
  echo 'Looked in metadata columns' "${PRIMER1_COLNUM}" 'and' "${PRIMER2_COLNUM}"
  echo 'Aborting script'
  exit
fi

# make primer array
read -a primers_arr <<< $( echo $PRIMER1 $PRIMER2 )

# Reverse complement primers
PRIMER1RC=$( echo ${PRIMER1} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )
PRIMER2RC=$( echo ${PRIMER2} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )


# make primer array
read -a primersRC_arr <<< $( echo $PRIMER1RC $PRIMER2RC )
echo 'PRIMER2RC=' "${PRIMER2RC}"
################################################################################
# Calculate the expected size of the region of interest, given the total size of fragments, and the length of primers and tags
################################################################################
#new param READ_LEN is length of read expected after primer removal
#echo "**** Expected Final Read Length for Cutadapt:", "${READ_LEN}"

################################################################################
# Find raw sequence files
################################################################################
# Look for any file with '.fastq' in the name in the parent directory
# note that this will include ANY file with fastq -- including QC reports!

LIBRARY_DIRECTORIES=$( find "$PARENT_DIR" -name '*.fastq*' -print0 | xargs -0 -n1 dirname | sort --unique )

# Count library directories and print the number found
# TODO if LIBRARY_DIRECTORIES is an array, its length is "${#LIBRARY_DIRECTORIES[@]}"
N_library_dir=$(echo $LIBRARY_DIRECTORIES | awk '{print NF}')
echo "${N_library_dir}"" library directories found:"

# Assign it to a variable for comparison
LIBS_FROM_DIRECTORIES=$(for i in $LIBRARY_DIRECTORIES; do echo "${i##*/}" ; done)

# Read library names from file or sequencing metadata
if [ "${READ_LIB_FROM_SEQUENCING_METADATA}" = "YES" ]; then
	LIB_COL=$(awk -F',' -v LIB_COL_NAME=$LIBRARY_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == LIB_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
	LIBS=$(awk -F',' -v LIBCOL=$LIB_COL 'NR>1 {print $LIBCOL}' $SEQUENCING_METADATA | uniq)
	N_libs=$(echo $LIBS | awk '{print NF}')
	echo "Library names read from sequencing metadata (""${N_libs}"") total"
	echo "${LIBS}"
else
	LIBS=$(tr '\n' ' ' < "${LIB_FILE}" )
	N_libs=$(echo $LIBS | awk '{print NF}')
	echo "Library names read from lib file (""${LIBS}"") total"
fi

# Check that library names are the same in the metadata and file system
if [ "$LIBS_FROM_DIRECTORIES" != "$LIBS" ]; then
	echo "Warning: Library directories and library names in metadata are NOT the same. Something will probably go wrong later..."
else
	echo "Library directories and library names in metadata are the same - great job."
fi


################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
# BOG data are on a network disk.  Copy file to mbonserver, process, then remove
## Do we need to do a library-level loop, I think all we want is all the .fastq files together?
################################################################################

for C_LIB in $LIBS; do
		echo "Library " "${PARENT_DIR}"/"${C_LIB}"
		#Make directory to hold sequence files
		LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${C_LIB}
		echo "${LIB_OUTPUT_DIR}"
		mkdir -p "${LIB_OUTPUT_DIR}"  #Make a directory named by library ID, which is unique

		raw_files=($( find "${PARENT_DIR}"/"${C_LIB}" -name '*.fastq*' ))
	 	for myfile in "${raw_files[@]}"; do
				echo $(date +%H:%M) 'copying ' "${myfile}" "${LIB_OUTPUT_DIR}"/"${myfile##*/}"
				cp "${myfile}"  "${LIB_OUTPUT_DIR}"/"${myfile##*/}"
		done

  	echo "find READS in " "${LIB_OUTPUT_DIR}"
 		# Identify the forward and reverse fastq files.
		READS=($(find "${LIB_OUTPUT_DIR}" -name '*.fastq*' | sort))
		READ1="${READS[0]}"
		READ2="${READS[1]}"

		##############################################################################
		# Atropos Remove Primers First from .fastq files
		##############################################################################
		# A change is that we're running this on each unmerged .fastq file prior to input into Dada2
		# following tutorial here: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#dada2


		################################################################################
		# PRIMER REMOVAL
		################################################################################
		# Atropos modification: PRIMERS are sometimes in reads and sometimes not depending on marker,
		# will not discard reads that do not have primers in them
		echo 'Remove Primers?: ' "${REMOVE_PRIMER}"
		if [ "${REMOVE_PRIMER}" = "YES" ]; then
				echo '[ "${REMOVE_PRIMER}" == "YES" ]'
				echo '**Forward read fastq=' "${READ1}"  #filename determined above
				echo '**Reverse read fastq=' "${READ2}"  #filename determined above
		    echo '**Results Files' "${TAG_DIR}"  #filename determined above
				#This removes primers sequentially instead of in the paired-end protocol which requires the presence of the first primer
				echo 'primer1' "${PRIMER1}"
				echo 'primer1RC' "${PRIMER1RC}"
				echo 'primer2' "${PRIMER2}"
				echo 'primer2RC' "${PRIMER2RC}"
				#Paired-end Processing
				#Atropos
				#atropos -a Rrev -g F -A Frev -G R
				#Change number of threads with -T
				atropos -a "${PRIMER2RC}" -g "${PRIMER1}" \
				    -A "${PRIMER1RC}" -G "${PRIMER2}"  \
						-o "${LIB_OUTPUT_DIR}"/"${C_LIB}"_R1_trimmed.fq.gz -p "${LIB_OUTPUT_DIR}"/"${C_LIB}"_R2_trimmed.fq.gz \
						--threads 15 \
						--times 10 \
						-e 0.20 \
						--minimum-length 20 \
						-pe1 "${READ1}" -pe2 "${READ2}"
				wait
			else
				#need to zip?
				cp "${READ1}" "${LIB_OUTPUT_DIR}"/"${C_LIB}"_R1_trimmed.fq
				cp "${READ2}" "${LIB_OUTPUT_DIR}"/"${C_LIB}"_R2_trimmed.fq
			fi

done


################################################################################
# END LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################


################################################################################
# Run Dada2 Over ALL Files
################################################################################

# RUN with Rscript Dada2_2019v1.0.r Dir_path
echo '#RUN with Rscript Dada2_2019v1.0.r Dir_path'
echo $(date +%Y%m%d_%H%M)
#e.g. : Rscript ./Dada2_2019v1.0.r /Users/kpitz/Projects/Dada2_Testing/CN18F_12S_2

echo "Run Dada2 R script: Rscript "${SCRIPT_DIR}"/Dada2_2019v1.0.r "${ANALYSIS_DIR}" "${truncF}" "${truncR}""

(cd "${ANALYSIS_DIR}" && Rscript "${SCRIPT_DIR}"/"Dada2_2019v1.0.r" "${ANALYSIS_DIR}" "${truncF}" "${truncR}")

################################################################################
# BLAST CLUSTERS with mbon-master spawning VMs
################################################################################
# This script uses GNU parallel to spawn several blast processes on different servers.
# After the processing has been completed, the worker will copy the file to the
# /tmp/Worker_Output directory on the mbon-master system, and the contents of this folder will
# be moved to /Blast_VM_output

BLAST_INPUT="${ANALYSIS_DIR}/ASVs.fa"

current_time=$(date +"%T")
echo "Blast started : $current_time"
echo "Number of ASVs:"
(grep '>' "$BLAST_INPUT" | wc -l)


blast_output="${STANDARD_PREFIX}_${DB_NAME}.xml"

echo ":::::::::::::::::::::: BLAST ::::::::::::::::::::::::::::::::"
echo $(date +%Y%m%d_%H%M) " BLASTing..."
echo " Database: " "${BLAST_DB}"
echo " Percent identity: " "${PERCENT_IDENTITY}"
echo " Word size: " "${WORD_SIZE}"
echo " E value: " "${EVALUE}"
echo " Maximum target matches: " "${MAXIMUM_MATCHES}"
echo " Output format: 5"
echo " Blast input: " "${BLAST_INPUT}"
echo " Blast output: " "${blast_output}"

# Temporary destination for output files
TempDestination="/tmp/Worker_Output"
#remove any files that may be present in this folder
rm -f ${TempDestination}/*
# Final destination for these output files after all the blasting is done
FinalDestination="${ANALYSIS_DIR}/Blast_VM_output"


# The number of cores for blast to use on each system
n_cores=8
# Number of virtual machines
NumWorkers=6

# Calculate the size of each of the file chunks
#cd $FileDir
# Get the number of records in the input file
NumRecords=`grep "^>" $BLAST_INPUT | wc -l`
echo "Number of Records:"
echo "$NumRecords"
# Get the size of each chunk depending on the number of workers
SizePerWorker=$(($NumRecords / $NumWorkers))
echo "Worker File size:"
echo "$SizePerWorker"

# Add 1 just to avoid rounding errors
SizePerWorker=$(($SizePerWorker + 1 ))
echo "Worker final file size:"
echo "$SizePerWorker"


cat $BLAST_INPUT | \
parallel -S mbon-worker1,mbon-worker2,mbon-worker3,mbon-worker4,mbon-worker5,mbon-worker6 \
--N $SizePerWorker \
--recstart '>' \
--pipe blastn-wrapper \
    -db "$BLAST_DB" \
    -num_threads "$n_cores" \
    -perc_identity "${PERCENT_IDENTITY}" \
    -word_size "${WORD_SIZE}" \
    -evalue "${EVALUE}" \
    -max_target_seqs "${MAXIMUM_MATCHES}" \
    -outfmt 5 \
    -out "${blast_output}" \
    -gapopen 5 \
    -gapextend 2 \
    -reward 2 \
    -penalty -3 \
    -query -



#output should be saved in analysis directory with different appended names for worker VMs
#Will have multiple xml files as output.

# Make sure the output directory exists
mkdir $FinalDestination > /dev/null 2>&1
# Move the output files from the temporary destination to the final destination
mv ${TempDestination}/* $FinalDestination

echo $(date +%Y%m%d_%H%M) " Blast finished."

# check for blast output directory and make sure it isn't empty
[ "$(ls -A "${FinalDestination}")" ] && echo "Blast Not Empty" || echo "Blast Empty"

#list blast files that were created
(ls -ltr ${FinalDestination})

#list number of unique query IDs that were blasted, make sure it matches number of ASVs
(grep '<Iteration_query-def>' ${FinalDestination}/*.xml | wc -l)


echo "Create Post_Blast_processing File"
################################################################################
# Create Post_Blast_processing File
################################################################################

if [ "$(ls -A "${FinalDestination}")" ]; then

	echo $(date +%H:%M) "${FinalDestination}" ',BLAST output found; proceeding to create Post_Blast file.'
	blast_fail="NO"
	POST_BLAST_FILE="${MY_SCRIPT_DIR}"/Post_BLAST_processing_v2.0_"$START_TIME".sh
	echo "POST_BLAST_FILE=" "${POST_BLAST_FILE}"
	echo SCRIPT_DIR="$SCRIPT_DIR" >> "$POST_BLAST_FILE"
	#echo RSCRIPTS="/home/mbonteam/dev" >> "$POST_BLAST_FILE"
	#
	echo "###################################################################################################" >> "$POST_BLAST_FILE"
	echo  STANDARD_PREFIX="$STANDARD_PREFIX" >> "$POST_BLAST_FILE"
	echo  DB_NAME="$DB_NAME" >> "$POST_BLAST_FILE"
	echo  BIOM_FILE_NAME="$BIOM_FILE_NAME" >> "$POST_BLAST_FILE"
	echo  ANALYSIS_DIR="$ANALYSIS_DIR" >> "$POST_BLAST_FILE"
	#echo  DIR="$ANALYSIS_DIR"/all_lib >> "$POST_BLAST_FILE"
	echo  START_TIME="$START_TIME" >> "$POST_BLAST_FILE"
	echo  param_file="$param_file" >> "$POST_BLAST_FILE"
	#echo  blast_output="$blast_output" >> "$POST_BLAST_FILE"
	echo  FinalDestination="$FinalDestination" >> "$POST_BLAST_FILE"
	echo  BLAST_INPUT="${BLAST_INPUT}" >> "$POST_BLAST_FILE"
	echo "###################################################################################################" >> "$POST_BLAST_FILE"
	cat /home/mbonteam/dev/Post_BLAST_processing_dada2_v2.0.sh  >> "$POST_BLAST_FILE"
else
	echo
	echo 'BLAST failed: the output file is empty or absent.'
	echo 'Files should be in:' "${FinalDestination}"
	echo
	blast_fail="YES"
fi # if blast_output


#<Iteration_query-ID> unique string per query in xml file
otu_cnt=$(wc -l "${ANALYSIS_DIR}"/ASVs_counts.tsv )
otu_fasta_cnt=$(grep -o  '>' "$BLAST_INPUT"  | wc -l)
blast_cnt=$(grep '<Iteration_query-def>' ${FinalDestination}/*.xml | wc -l)


echo $(date +%H:%M) "Compressing output fasta files..."
find "${ANALYSIS_DIR}" -type f -name '*.fa' -exec ${ZIPPER} "{}" \;
echo $(date +%H:%M) "Removing original and Dada2 filtered fastq files..."
rm -r "${ANALYSIS_DIR}"/filtered  #get rid of the Dada2 filtered fastq files
rm "${ANALYSIS_DIR}"/*/*.fastq.gz   #catch the original fastq files in library folders


#need xml files unzipped for Post Blast processing.
#find "${ANALYSIS_DIR}" -type f -name '*.xml' -exec ${ZIPPER} "{}" \;
echo $(date +%H:%M) "Cleanup performed."


if [ "$blast_fail" = "NO" ]; then
	echo -e '\n'$(date +%H:%M)"\t Pipeline finished! Here's a Terry Pratchett quote:" '\n'
	echo -e '\tWhy do you go away? So that you can come back. So that you can see the place '
	echo -e '\tyou came from with new eyes and extra colors. And the people there see you differently, too.'
	echo -e '\tComing back to where you started is not the same as never leaving.'
	echo -e '\t― Terry Pratchett, A Hat Full of Sky''\n'

	echo '******************** TAXONOMIC ASSIGNMENT **********************'
	echo '*                                                              *'
	echo '*   Run Post_BLAST_processing file from your script directory  *'
	echo '*        bash  Post_BLAST_processing_v2.0_'${START_TIME}'.sh        *'
	echo '*                                                              *'
	echo '*********************** MEGAN & BIOM ***************************'

fi


echo ''
echo ''
echo 'OTU_table count=' "$otu_cnt"
echo 'OTU_fasta count=' "$otu_fasta_cnt"
echo 'BLAST count=' "$blast_cnt"

FINISH_TIME=$(date +%Y%m%d_%H%M)

if [ "$NOTIFY_EMAIL" = "YES" ]; then
	echo 'Dada2 Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME 'Data here:' "${ANALYSIS_DIR}" | mail -s "banzai is finished" "${EMAIL_ADDRESS}"
else
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME
fi
echo 'Data are in ' "${ANALYSIS_DIR}"

#code to create a simplified logfile - need to revise with current output.
# remove some of the repetitive lines, read percents from kmers, fastq etc.  makes it difficult to read logfile.
# grep -v flag means output everything except pattern, -f means pattern is in text file, i.e. grep_patterns.txt
#mv "${ANALYSIS_DIR}"/logfile.txt "${ANALYSIS_DIR}"/logfile_orig.txt
#grep -v -f /home/mbonteam/MBARI/reiko/scripts/grep_patterns.txt "${ANALYSIS_DIR}"/logfile_orig.txt > "${ANALYSIS_DIR}"/logfile.txt
# Quick scan of logfile2.txt for errors, error line and 3 lines before, -B 3
grep -n -B 3 "error" "${ANALYSIS_DIR}"/logfile.txt >"${ANALYSIS_DIR}"/errors.txt
