#!/usr/bin/env python

'''
Script written by Katie Pitz 12/13/19
Takes Dada2 output and MEGAN6 output and formats into csv files:
OTU_table
Taxa_table

Call with
python Make_otu_taxa_table.py ANALYSIS_DIR

Updated 072523 by kpitz to better parse MEGAN results
'''
import sys
import pandas as pd
import numpy as np
import glob

#Get analysis directory from Post_Blast_Processing.sh script
directory = sys.argv[1] +'/'

#Get File paths
meta_file = glob.glob(directory+'*_analysis_metadata.csv')[0]
print('Metadata:', meta_file)
otu_table_file = glob.glob(directory+'ASVs_counts.tsv')[0]
print('ASV table:', otu_table_file)
taxa_table_file = glob.glob(directory+'tpath.txt')[0]
print('taxa table:', taxa_table_file)

#Import Metadata
df = pd.read_csv(meta_file)
df.set_index('library', inplace=True)
meta = df.copy()
biom_meta = meta.fillna('Nan')
biom_meta.reset_index(inplace=True)
biom_meta.rename(columns={'sample_name':'#sample_name'}, inplace=True) #Biom conversion needs # in index header
biom_meta.set_index('#sample_name', inplace=True)

#export biom meta with no empty fields and sample_name as index:
biom_meta.to_csv(directory + 'Biom_metadata.tsv', sep='\t')

#OTU table
df = pd.read_csv(otu_table_file, sep='\t')
df.rename(columns={'Unnamed: 0':'ASV'}, inplace=True)
df.set_index('ASV', inplace=True)
otu_table = df.copy()

#Taxa table
#levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
#updated 072523 kpitz
levels = ['Domain', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
df = pd.read_csv(taxa_table_file, sep="\t", header=None, names =['ASV', 'taxonomy'])

# extract all levels present in the taxonomy:
df.set_index('ASV', inplace=True)

df['Domain'] = df['taxonomy'].str.extract(r'\[D\] ([^;]*);.*')
df['Kingdom'] = df['taxonomy'].str.extract(r'\[K\] ([^;]*);.*')
df['Phylum'] = df['taxonomy'].str.extract(r'\[P\] ([^;]*);.*')
df['Class'] = df['taxonomy'].str.extract(r'\[C\] ([^;]*);.*')
df['Order'] = df['taxonomy'].str.extract(r'\[O\] ([^;]*);.*')
df['Family'] = df['taxonomy'].str.extract(r'\[F\] ([^;]*);.*')
df['Genus'] = df['taxonomy'].str.extract(r'\[G\] ([^;]*);.*')
df['Species'] = df['taxonomy'].str.extract(r'\[S\] ([^;]*);.*')

df = df.sort_values(levels)

# MEGAN already exports 'unknown' term but is inconsistent for some terms, like 'Kingdom'
# if no taxonomy:
for level in levels:
    df.loc[df['taxonomy'].isna(), level] = 'unassigned'

df.drop('taxonomy', axis=1, inplace=True)
taxa_tab = df.copy()

#Export Files
df = pd.concat([meta[['sample_name']], otu_table.T], axis=1, sort=False)
#rename columns as sample_names
df.set_index('sample_name', inplace = True)
#Join with taxonomy
samples = list(df.T)
print('Number of Samples in OTU table:', len(samples))
df = pd.concat([taxa_tab, df.T], axis=1, sort=False)
df.fillna('no_hit', inplace=True)
#sort by ASV number - relates to abundance of ASV
df['OTU_num'] = df.index.str.split('_').str[-1].astype(int)
df=df.sort_values('OTU_num')
df.drop('OTU_num', inplace=True, axis=1)
#Make index names
df.index.name='ASV'
#export 'OTU_taxa_table_all.csv'
df.to_csv(directory+'ASV_taxa_table_all.csv')
#export taxa table
taxa_tab_all = df[levels].copy()
taxa_tab_all.index.name='#ASV'  #Biom conversion needs # in index header
taxa_tab_all.to_csv(directory+'Taxa_table.tsv', sep='\t')
#export otu table
otu_tab_all = df[samples].copy()
otu_tab_all.to_csv(directory+'ASV_table.tsv', sep='\t')
#Print some basic stats of total reads and ASVs
df['Total_Reads'] = df[samples].sum(axis=1)
df['Total_ASVs']=1
df=df.groupby(['Kingdom','Phylum']).sum()
print('XXXXXX Summary of Total Reads and ASVs by Kingdom and Phylum')
print(df[['Total_Reads','Total_ASVs']])
print('XXXXXX End Summary')
