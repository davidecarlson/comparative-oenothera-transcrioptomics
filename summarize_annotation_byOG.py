#!/usr/bin/env python

import glob
import os
import re
import sys
from itertools import chain
import pandas as pd
from collections import Counter
import csv
import multiprocessing as mp
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag


# path to trinotate functional annotation tables

trinotate_path = '/datahome/oenothera/transcripts/assemblies/snakemake-trinotate/results/reports/'

# read in file with GO terms

godag = GODag('go-basic.obo',optional_attrs={'relationship'},load_obsolete=True)
#optional_relationships = set()

# lists of GO umbrella terms

MF_umb = ['GO:0005488', 'GO:0003824', 'GO:0005215', 'GO:0045182', 
'GO:0140110', 'GO:0005198', 'GO:0060089', 'GO:0016209',
'GO:0098772', 'GO:0060090', 'GO:0045735', 'GO:0038023', 'GO:0016247',
'GO:0140104', 'GO:0038024','GO:0060090']

BP_umb = ['GO:0065007', 'GO:0000003', 'GO:0040011', 'GO:0009987', 'GO:0032502', 
'GO:0051179', 'GO:0008152', 'GO:0050896', 'GO:0023052', 'GO:0016043', 'GO:0001906', 'GO:0022610']

CC_umb = ['GO:0016020', 'GO:0000785', 'GO:0009295', 'GO:0030054', 'GO:0032991', 'GO:0005737', 
'GO:0031974','GO:0110165','GO:0031975', 'GO:0019012', 'GO:0005576', 'GO:0043226', 'GO:0045202', 'GO:0099080']


#OG_file = 'sample_OGs.tsv'
OG_file = 'Orthogroups_noheader.tsv'

file = open(OG_file, 'r')
OG = file.read().replace(' ', '')


newOG = OG.replace(' ', '').replace(',', '\t').strip()

OG_list = [line.split('\t') for line in newOG.split('\n') if line]

num_OGs = len(OG_list)
print(num_OGs)


# create function that will take an OG and return all GO annotations for all the proteins in that OG
def get_annot_results(OG):

	OG_ID = OG[0]


	print("Finding annotation results for orthogroup " + OG_ID)

	proteins = OG[1:]
	proteins =  list(filter(None, proteins))


	seen_CCs = set()
	seen_MF = set()
	seen_BP = set()

	CC_list = []
	CC_dict = {}
	MF_list = []
	BP_list = []

	for protein in proteins:

		name = re.split(r'(-[A-Z]{4}_)', protein)[2]
		sample_ID = protein.split('-')[1].split('_')[0]

		#print("Looking up functional annotations")
		with open(trinotate_path + sample_ID + ".trinotate_annotation_report.tsv") as anno_file:
			for line in anno_file:
				if line.split('\t')[4]==name and line.split('\t')[13] != '.':

					#print(line)
					blastx_annot = line.split('\t')[13]
					annot_list = blastx_annot.split('`')
					#print(annot_list)
					CC = [i.split('^')[0] for i in annot_list if i.split('^')[1].startswith('cellular_component')]
					MF = [i.split('^')[0] for i in annot_list if i.split('^')[1].startswith('molecular_function')]
					BP = [i.split('^')[0] for i in annot_list if i.split('^')[1].startswith('biological_process')]

					# check if annotation has already been added to the list for that protein from another transcript to reduce redundnancy
					for item in CC:
						if item not in seen_CCs:
							seen_CCs.add(item)
							CC_list.append(CC)

					for item in MF:
						if item not in seen_MF:
							seen_MF.add(item)
							MF_list.append(MF)
					for item in BP:
						if item not in seen_BP:
							seen_BP.add(item)
							BP_list.append(BP)

	# flatten the list of lists and remove duplicates
	CC_list = list(chain(*CC_list))
	MF_list = list(chain(*MF_list))
	BP_list = list(chain(*BP_list))

	# remove duplicate entries from the list
	CC_list = list(set(CC_list))
	MF_list = list(set(MF_list))
	BP_list = list(set(BP_list))

	return(CC_list,MF_list, BP_list)

# parallelize finding all the annotations over each of the OGs

with mp.Pool(processes=40) as pool:
	results = pool.starmap_async(get_annot_results, zip(OG_list))
	myresults = results.get()

pool.close()
pool.join() # Wait for all child processes to close.

#print(myresults)

# split results by the 3 GO categories and flatten the lists
myCC = list(chain(*[result[0] for result in myresults if result[0]]))
myMF = list(chain(*[result[1] for result in myresults if result[1]]))
myBP = list(chain(*[result[2] for result in myresults if result[2]]))


# count how many times each GO annotation appears and get the proportion of OGs with annotation

def count_annot(annot_list, length):
	annot_counter = Counter(annot_list)
	results = []
	for key, value in annot_counter.items():
		count = value
		prop = float(count) / length
		annot = key
		combined = [annot, count, prop]
		results.append(combined)
		#print(OG, count, prop)
	return(results)

CC_results = count_annot(myCC, num_OGs)
MF_results = count_annot(myMF, num_OGs)
BP_results = count_annot(myBP, num_OGs)


def lookup_terms(annot_list, umb_list):
	# for each GO ID, return the annotation name and any ancestors if they're in the list of umbrella terms
	umbrella_results = []
	for term in annot_list:
		#get the name associated with the ID
		name = godag[term[0]].name
		# find the ancestors of the GO term
		optional_relationships = {'regulates', 'negatively_regulates', 'positively_regulates'}
		subset = GoSubDag([term[0]], godag, relationships=optional_relationships, prt=None)
		ancs = list(subset.rcntobj.go2ancestors[term[0]])
		# check if any of the ancestors are in the umbrella term list
		umbrellas =  [item for item in ancs if item in umb_list]
		#print(umbrellas)
		umbrella_names = ','.join([godag[item].name for item in umbrellas])
		if not umbrella_names:
			if term[0] in umb_list:
				umbrella_names = godag[term[0]].name
			else:
				umbrella_names = 'none'
		results = [term[0], name,umbrella_names, term[1], term[2]]
		umbrella_results.append(results)
	return(umbrella_results)


CC_umbrellas = lookup_terms(CC_results,CC_umb)
MF_umbrellas = lookup_terms(MF_results,MF_umb)
BP_umbrellas = lookup_terms(BP_results,BP_umb)

# make dataframes and write results to files

df_CC_umbrellas = pd.DataFrame(CC_umbrellas, columns = ['GO Annotation ID', 'GO Annotation Name', 'Umbrella Annotation Names', 'Count of OGs with Annotation', 'Proportion of OGs with Annotation'])
df_CC_umbrellas.sort_values(by='Count of OGs with Annotation', ascending=False, inplace=True)
print(df_CC_umbrellas)
df_CC_umbrellas.to_csv('cellular_component_umbrella_results.csv', index = False)

df_MF_umbrellas = pd.DataFrame(MF_umbrellas, columns = ['GO Annotation ID', 'GO Annotation Name', 'Umbrella Annotation Names', 'Count of OGs with Annotation', 'Proportion of OGs with Annotation'])
df_MF_umbrellas.sort_values(by='Count of OGs with Annotation', ascending=False, inplace=True)
print(df_MF_umbrellas)
df_MF_umbrellas.to_csv('molecular_function_umbrella_results.csv', index = False)

df_BP_umbrellas = pd.DataFrame(BP_umbrellas, columns = ['GO Annotation ID', 'GO Annotation Name', 'Umbrella Annotation Names', 'Count of OGs with Annotation', 'Proportion of OGs with Annotation'])
df_BP_umbrellas.sort_values(by='Count of OGs with Annotation', ascending=False, inplace=True)
print(df_BP_umbrellas)
df_BP_umbrellas.to_csv('biological_process_umbrella_results.csv', index = False)
