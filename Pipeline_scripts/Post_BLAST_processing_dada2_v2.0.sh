#Begin Post_Blast_processing script v2.0
#Katie Pitz 121619
#Edited by Markus Min 091820 to include decontamination and standard outputs R script

OTU_table="${ANALYSIS_DIR}"/ASVs_counts.tsv
START_T=$(date +%Y%m%d_%H%M)
PDIR="${ANALYSIS_DIR}"/Post_Blast_"${START_T}"
mkdir "${PDIR}"

# Write a log file of output from this script (everything that prints to terminal)
LOGFILE="${PDIR}"/Post_BLASTlogfile.txt
exec > >(tee "${LOGFILE}") 2>&1
echo "Running script from: "
echo "$BASH_SOURCE"
echo "Parameters read from: "
echo "${param_file}"
source "${param_file}"
echo "metadata read from: "
echo "${SEQUENCING_METADATA}"
echo 'Post Blast Directory:'
echo "${PDIR}"
echo 'Script Directory:'
echo "${SCRIPT_DIR}"
echo '_____________________________________________________'
cp $BASH_SOURCE "${PDIR}"/Post_BLAST_analysis_script.txt

#Blast output is now several xml files within FinalDestination="${ANALYSIS_DIR}/Blast_VM_output"
#For each file in the folder, need to run through megan and get a rma6 file
#And then convert each rma6 file into a tab deliminated file, and cat these together
echo blast output is in "${ANALYSIS_DIR}/Blast_VM_output"
#For some reason, megan puts the output files one directory above the one you list.
#Seems like a bug in the program

otu_cnt=$(wc -l "${ANALYSIS_DIR}"/ASVs_counts.tsv )
#otu_fasta_cnt=$(grep -o  '>' "$BLAST_INPUT"  | wc -l)
blast_cnt=$(grep '<Iteration_query-def>' "${FinalDestination}"/*.xml | wc -l)


echo $(date +%H:%M) "${FinalDestination}" ',BLAST output; proceeding to MEGAN6.  Taxonomic ranks exported to files.'

echo /usr/local/megan/tools/blast2rma --in "${ANALYSIS_DIR}/Blast_VM_output/*.xml" --format BlastXML --blastMode BlastN --out "${ANALYSIS_DIR}/Blast_VM_output"/ --minScore "${MINIMUM_SCORE}" --maxExpected "${MAX_EXPECTED}" --topPercent "${TOP_PERCENT}" --minSupport "${MINIMUM_SUPPORT}" --lcaAlgorithm naive --lcaCoveragePercent "${LCA_PERCENT}"
/usr/local/megan/tools/blast2rma --in "${ANALYSIS_DIR}"/Blast_VM_output/*.xml --format BlastXML --blastMode BlastN --out "${ANALYSIS_DIR}/Blast_VM_output"/ --minScore "${MINIMUM_SCORE}" --maxExpected "${MAX_EXPECTED}" --topPercent "${TOP_PERCENT}" --minSupport "${MINIMUM_SUPPORT}" --lcaAlgorithm naive --lcaCoveragePercent "${LCA_PERCENT}"

#Now there should be a list of .rma6 files in "${ANALYSIS_DIR}/
echo 'Now for each Megan6 file, extract the taxonomy:'
for file in "${ANALYSIS_DIR}"/*.rma6
do
  echo "$file"
	echo /usr/local/megan/tools/rma2info --in "${file}"  --read2class Taxonomy --paths --majorRanksOnly > "${file}"_tpath.txt
	/usr/local/megan/tools/rma2info --in "${file}"  --read2class Taxonomy --paths --majorRanksOnly > "${file}"_tpath.txt
done


#join files together; get all taxonomy results
cat "${ANALYSIS_DIR}"/*_tpath.txt > "${ANALYSIS_DIR}"/tpath.txt

################################################################################
# FORMAT TAXONOMY RESULTS, OTU Tables
################################################################################
#call python script with location of Analysis folder
python3 "${SCRIPT_DIR}"/Make_otu_taxa_table.py "${ANALYSIS_DIR}"

#ASV_table.csv
#Taxa_table.csv
#ASV_taxa_table_all.csv


################################################################################
# Filter Taxonomy by Blast hit
################################################################################
# ::::::::::::::::::: Filter XML File for PhiX and apply genus, species limits ::::::::::::::::::::::::
#Filter out OTUs with hits to PhiX - remove from complete OTU table (will be unannotated)
#Filter annotations at the genus and species level to different bitscore, %ID limits
# search by genbank sequence name
#Call xml file -> feed into python script

echo ::: Filter XML File for PhiX and apply genus, species limits :::

otu_all="${ANALYSIS_DIR}"/ASV_taxa_table_all.csv

echo Blast XML files in : "${FinalDestination}"
echo OTU_table_taxa_all.txt :  "${otu_all}"
echo python3 "${SCRIPT_DIR}"/XML_Filter_dada2.py "${bitscore_sp}" "${per_ID_sp}" "${bitscore_gn}" "${per_ID_gn}" "${FinalDestination}" "${otu_all}"
python3 "${SCRIPT_DIR}"/XML_Filter_dada2.py "${bitscore_sp}" "${per_ID_sp}" "${bitscore_gn}" "${per_ID_gn}" "${FinalDestination}" "${otu_all}"


################################################################################
# CREATE BIOM FILE
################################################################################
echo :::::::::::::::::::::::: Create JSON BIOM file ::::::::::::::::::::::::
#Have to add tables in succession (otu table, then taxa table, then metadata)
biom convert -i "${ANALYSIS_DIR}"/ASV_table.tsv -o "${ANALYSIS_DIR}"/"${BIOM_FILE_NAME}".biom --table-type="OTU table" --to-json
biom add-metadata -i "${ANALYSIS_DIR}"/"${BIOM_FILE_NAME}".biom -o "${ANALYSIS_DIR}"/"${BIOM_FILE_NAME}"_obs.biom --observation-metadata-fp "${ANALYSIS_DIR}"/Taxa_table.tsv --output-as-json
biom add-metadata -i "${ANALYSIS_DIR}"/"${BIOM_FILE_NAME}"_obs.biom -o "${ANALYSIS_DIR}"/"${BIOM_FILE_NAME}"_obs_meta.biom --sample-metadata-fp "${ANALYSIS_DIR}"/Biom_metadata.tsv --output-as-json

################################################################################

#gzip "${BLAST_INPUT}"
################################################################################
echo "Create a PEAR report with the names of the read files and the assembled, discarded, and not assembled read numbers.  The report also includes the PEAR parameters"
echo 'python3 '"${SCRIPT_DIR}"'/mk_PEAR_report.py ' "${ANALYSIS_DIR}"
python3 "${SCRIPT_DIR}"/mk_PEAR_report.py "${ANALYSIS_DIR}"

FINISH_TIME=$(date +%Y%m%d_%H%M)

################################################################################

# RUN R SCRIPT FOR DECONTAMINATION + STANDARD OUTPUTS
################################################################################
INPUT_ASV_TABLE="${ANALYSIS_DIR}"/Filtered_ASV_taxa_table_all.csv
METADATA=$(find "${ANALYSIS_DIR}" -name '*_analysis_metadata.csv')

echo $METADATA
echo $INPUT_ASV_TABLE

Rscript "${SCRIPT_DIR}"/decon_std_outputs_v1.0.R $INPUT_ASV_TABLE $METADATA "${ANALYSIS_DIR}"


if [ "$NOTIFY_EMAIL" = "YES" ]; then
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME | mail -s "banzai is finished" "${EMAIL_ADDRESS}"
else
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME
fi
echo 'Data are in ' "${PDIR}"


echo ''
echo ''
echo 'OTU_table count=' "$otu_cnt"
#echo 'OTU_fasta count=' "$otu_fasta_cnt"
echo 'BLAST count=' "$blast_cnt"
#echo 'Megan count=' "$megan_cnt"
#


echo ''
echo ''
echo -e '\n'$(date +%H:%M)'\tBiom finished! Time for a G&T...\n'
echo
echo -e '\t~~~~ ><(((°> Gin & Tonic <°)))>< ~~~~'
echo -e '\t2 oz\t Gin'
echo -e '\t5 oz\tChilled Tonic Water'
echo -e '\t1\tSlice of Lime'
echo -e '\t2\tCubes of Ice'
echo -e '\t1\tChilled Martini Glass'
echo -e '\tShake ice and gin until well chilled'
echo -e '\tStrain into chilled Martini glass'
echo -e '\tAdd tonic and lime. '
echo -e '\tPut your feet up and relax, life is good.' '\xf0\x9f\x8d\xb9\x0a''\n'




( set -o posix ; set ) > "${PDIR}"/post_blast_variables.txt

#move files to Post_BLAST DIR
mv   "${ANALYSIS_DIR}"/*rma6*  "${PDIR}"/
mv   "${ANALYSIS_DIR}"/*tpath*  "${PDIR}"/


# Quick scan of Post_BLASTlogfile.txt for errors
#source "${SCRIPT_DIR}"/post_blast_errors.sh  "${PDIR}" "${BASH_SOURCE}"
