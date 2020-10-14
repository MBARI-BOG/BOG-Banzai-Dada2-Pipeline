#!/usr/bin/env python

'''
Written by Katie Pitz
Updated 121119 to work on multiple xml files produced by 'mbon-master' VM system KP

updated 04/10/18
XML_Filter.py : Go through BLAST xml results file and pull out OTUs which have a hit to
PhiX. Also impose limits on bitscore and Percent Identity for OTU assignment at the
species or genus level.

From New_Post_Blast.sh:
#BLAST parameters for species assignment:
#18S bitscore ~ 200 ; COI bitscore ~ 400 (COI has a longer alignment length)
bitscore_sp = 500
per_ID_sp = 97
#BLAST parameters for genus assignment:
#18S bitscore ~ 150 ; COI bitscore ~ 400 (COI has a longer alignment length)
bitscore_gn = 400
per_ID_gn = 95

call with:
python XML_Filter.py bitscore_sp per_ID_sp bitscore_gn per_ID_gn XML_File.xml otu_taxa_table_all.txt

product is a new OTU table with taxa annotations bumped up and no PhiX-hitting OTUs
blast stats csv file
'''

import glob
import sys
from Bio.Blast import NCBIXML
import pandas as pd


#Parameters in sys.argv
#Species limits
bitscore_sp = int(sys.argv[1])
per_ID_sp = int(sys.argv[2])
print('Species bitscore limit: ', bitscore_sp, ', Percent Identity Limit: ', per_ID_sp)
#Genus limits
bitscore_gn = int(sys.argv[3])
per_ID_gn = int(sys.argv[4])
print('Genus bitscore limit: ', bitscore_gn, ', Percent Identity Limit: ', per_ID_gn)

# XML File locations
directory = sys.argv[5]
files = glob.glob(directory+'/*.xml')
#print(files)

#Import OTU tables
otu_file = sys.argv[6]
df = pd.read_csv(otu_file, sep=',')

df.set_index('ASV', inplace=True)
levels = ['Kingdom','Phylum','Class','Order','Family','Genus','Species']

#OTU Table
otu_tab = df.drop(levels, axis=1)

#Taxa Table
taxa_tab = df[levels]

#Create dictionary of ASV ID to taxonomy
taxa_dict = taxa_tab.T.to_dict('list')

# Parse through XML files
xml_dfs = []  #will store results from each xml file
for filename in files:
    #parse filename and directory location
    print('Blast XML file: ', filename)
    directory = filename.split('/')[:-2]
    string_s= '/'
    directory=string_s.join(directory)
    directory = directory +'/'
    print('directory', directory)
    genus=[]
    species=[]
    query=[]
    hit_dict ={}  #Dictionary to store hits
    result_handle =open(filename)
    #Phix names to filter out - added for diversity in sequencing process
    taxon = ['PhiX', 'phiX']
    #parse blast records
    blast_records= NCBIXML.parse(result_handle)
    #Locations to store PhiX hits
    taxon_A =[]
    taxon_B =[]
    query_A =[]
    query_B =[]
    for blast_record in blast_records:
        query = blast_record.query
        genus = str(taxa_dict[query][-2]) #look up genus assignment
        species = str(taxa_dict[query][-1]) #look up species assignment
        hit_ID=0 #hit counter
        for alignment in blast_record.alignments:
            hit_ID+=1
            key = (query, hit_ID)
            hsp = alignment.hsps[0] #only look at top hsp per alignment to genbank sequence
            #set limits on evalue and bitscore and % Identity
            per_iden = hsp.identities/float(hsp.align_length)
            per_iden = (per_iden *100)
            per_iden = round(per_iden, 2)
            evalue = hsp.expect
            bitscore = int(hsp.score)
            #alignment length
            align_len = hsp.align_length
            #Check if MEGAN assigned ASV at species level;look for species in hit ID if yes
            if species in ['s_', 'no_hit', 'unassigned']:
                species_h = 'species_unassigned'
            else:
                if species in alignment.hit_def:
                    species_h = species
                else:
                    species_h ='Not_Found'
            #Check if MEGAN assigned ASV at genus level;look for species in hit ID if yes
            if genus in ['g_', 'no_hit', 'unassigned']:
                genus_h = 'genus_unassigned'
            else:
                if genus in alignment.hit_def:
                    genus_h = genus
                else:
                    genus_h ='Not_Found'
            value = (evalue, bitscore, per_iden, genus_h, species_h, alignment.hit_def, align_len)
            hit_dict[key]=value
            genus_h='' #reset variables
            species_h=''#reset variables
            #Look for PhiX hits
            if taxon[0] in alignment.title:
                #print (alignment.title)
                taxon_A.append(blast_record.query)
                query_A.append(alignment.title)
            if taxon[1] in alignment.title:
                #print (alignment.title)
                taxon_B.append(blast_record.query)
                query_B.append(alignment.title)

    print ('Done parsing ', filename)

    #Save Blast Results to file
    df = pd.DataFrame(hit_dict)
    df=df.T
    df.columns = ['eval','bitscore', '%ID', 'genus', 'species', 'hit_def', 'align_len']
    df.reset_index(inplace=True)
    df = df.rename(columns={'level_0': 'OTU', 'level_1': 'Hit_number'})
    #Keep the best hit of each category of hit
    df = df.drop_duplicates(subset=['genus', 'species', 'OTU'], keep='first')
    #Sort by OTU number
    df['OTU_Number']= df['OTU'].str.split('_').str[1]
    df['OTU_Number'] = df['OTU_Number'].astype(int)
    df.sort_values(['OTU_Number','Hit_number'], inplace=True)
    df.drop('OTU_Number', axis=1, inplace=True)
    df.set_index(['OTU', 'Hit_number'], inplace = True)
    outfile = filename.replace('.xml','_TopBLASThits.csv')
    print(outfile)
    df.to_csv(outfile)
    xml_dfs.append(df)

#Look at all XML Files together

df= pd.concat(xml_dfs, axis=0)

# Begin filtering Taxonomic Annotations by Species and Genus limits
df[['eval', 'bitscore', '%ID', 'align_len']] = df[['eval', 'bitscore', '%ID', 'align_len']].apply(pd.to_numeric)
#Get top hits for taxa annotation assigned by MEGAN
df = df.reset_index().set_index('OTU')
df = df.join(taxa_tab, how='left')   #Join blast results with taxa table; anything without a blast hit won't join
df.reset_index(inplace=True)
df=df.sort_values(['index', 'Hit_number'])
df.set_index(['index', 'Hit_number'], inplace=True)

#Need to limit by top hit that includes genus and species, or just genus if that's all that's available
#Or if neither of them are available, set g_ and s_
#need to sort df by these levels; higher number, less specific assignment
df['G_lev']= 1 # 1 = assigned at that level; start out with everything =1
df['S_lev']= 1 # 1 = assigned at that level; start out with everything =1
# Now find unassigned ASVs:
df.loc[df['genus'] == 'genus_unassigned' , 'G_lev'] = 3
df.loc[df['genus'] == 'Not_Found' , 'G_lev'] = 2
df.loc[df['species'] == 'species_unassigned' , 'S_lev'] = 3
df.loc[df['species'] == 'Not_Found' , 'S_lev'] = 2

df = df.sort_values(['S_lev', 'G_lev'])
#Now just want to take the top hit for each OTU.
df=df.reset_index()
df = df.drop_duplicates('index', keep='first')
df.drop(['G_lev', 'S_lev'], axis=1, inplace=True)

#check no duplicates, one entry per OTU
dups = df.reset_index().duplicated(subset=['index'])
if dups.unique() != [False]:
    print('Error!! Duplicate OTUs in Table. Something went wrong.')


#Limit by percent ID
df.loc[df['%ID'] < per_ID_gn , 'Genus'] = "g_"
df.loc[df['%ID'] < per_ID_sp , 'Species'] = "s_"
#Limit by bitscore
df.loc[df['bitscore'] < bitscore_gn , 'Genus'] = "g_"
df.loc[df['bitscore'] < bitscore_sp , 'Species'] = "s_"

#sort by ASV number and hit number
df.reset_index(inplace=True)
df['OTU_Number']= df['index'].str.split('_').str[1]
df['OTU_Number'] = df['OTU_Number'].astype(int)
df.sort_values(['OTU_Number','Hit_number'], inplace=True)
df.set_index(['index','Hit_number'], inplace=True)

#Create New Taxa Table from selection
taxa_lim = df[levels]  #limit columns to taxonomy
taxa_lim = taxa_lim.reset_index().drop('Hit_number', axis=1)
taxa_lim = taxa_lim.rename(columns = {'index':'ASV'})
taxa_lim.set_index('ASV', inplace=True)

outname = directory + 'Filtered_ASV_taxa_table.csv'
print(outname)
taxa_lim.to_csv(outname)

#Get average %ID of genus and species corrected ASVs:
print('XXXX  Check Stats of hits that have been filtered: XXXX')
df_g = df.loc[df['Genus']=='g_']
print('Mean %ID of ASVs with g_:',df_g['%ID'].mean())
print('Mean bitscore of ASVs with g_:',df_g['bitscore'].mean())
df_s = df.loc[df['Species']=='s_']
df_s = df_s.loc[df_s['Genus']!='g_']
print('Mean %ID of ASVs with s_:',df_s['%ID'].mean())
print('Mean bitscore of ASVs with s_:',df_s['bitscore'].mean())

df.drop('OTU_Number', axis=1, inplace=True)
outname = directory + 'ASV_tophit_Ftaxonomy.csv'
print(outname)
df.to_csv(outname)


#Remove PhiX OTUs
list_DUP = (taxon_A, taxon_B)
#Number of OTUs without duplicates
for i in range(len(list_DUP)):
    y=list(set(list_DUP[i]))   #Remove duplicate entries from each list
    print(taxon[i], ' number of OTUs: ', len(y))    #Get number of OTUs which hit each taxon

#Get total number of otus and reads with PhiX and phiX; list of unique ASVs with PhiX hits
taxons = taxon_A + taxon_B
taxons= list(set(taxons))
print('Total Number OTUs with PhiX hits:',len(taxons))   #Number of OTUS

if len(taxons)== 0 :
    print('No Phix OTUs removed')
else:
    #remove these ASVs from the OTU table
    df = otu_tab.copy()
    df.reset_index(inplace=True)
    df=df.loc[df['ASV'].isin(taxons)==False]
    otu_tab = df.copy()

#ASV taxa table taken from Blast XML files will only contain ASVs with blast hits
#Ideally want to assign taxonomy from old table for those hits that weren't changed by this filtering
#They should all be 'no_hit'

#Join  new taxonomy and old taxonomy table together; check 'no_hit' for all missing taxa
df= taxa_lim.copy()
df['New']=1
df=pd.concat([df, taxa_tab], axis=1, keys=['new','old'],sort=False)
df=df.sort_values(('new', 'New'))
df=df.loc[df[('new', 'New')]!=1]
print('Unique Kingdom IDs from taxa without a Blast hit: ')
print("Should just be 'no_hit':")
print(df[('old', 'Kingdom')].unique())

#Since they're all no_hit, assign any missing values as 'no_hit'
df=pd.concat([taxa_lim, otu_tab], axis=1,sort=False)
df=df.fillna('no_hit')
df.index.name ='ASV'

#export 'OTU_taxa_table_all.csv'
outname = directory + 'Filtered_ASV_taxa_table_all.csv'
print(outname)
df.to_csv(outname)

print('Finished Filtering Taxonomy by BLAST hits')
