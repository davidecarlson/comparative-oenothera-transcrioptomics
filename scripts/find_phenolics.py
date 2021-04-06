#!/usr/bin/env python

import glob
from itertools import dropwhile
import re
import multiprocessing as mp
from itertools import chain
import pandas as pd

# search through each of the Trinotate transcript annotation reports and find transcripts with blastx matches to phenolic proteins

#list of phenolics protein names

phenolics = ['Flavonol synthase/flavanone 3-hydroxylase','Trans-cinnamate 4-monooxygenase', 'Phenylalanine ammonia-lyase',
                        '4-coumarate--CoA ligase','Chalcone synthase','Chalcone--flavonone isomerase',
			'Flavanone 3-dioxygenase','Naringenin,2-oxoglutarate 3-dioxygenase',"Flavonoid 3'-monooxygenase",
			"Flavonoid 3',5'-hydroxylase"]

#print(phenolics)

# create phenolics regex pattern for matching later

pattern = re.compile(r'\b(?:%s)\b' % '|'.join(phenolics))

# make dictionary of ID:taxon mappings

species_map = { 'GVCB': 'affinis-GVCB',
                'OEAR': 'argillicola-OEAR',
                'EQYT': 'berlandieri-EQYT',
                'EXGW': 'biennis-EXGW',
                'KBRW': 'clelandii-KBRW',
                'AJBK': 'elata-AJBK', 
                'OEEE': 'elata-OEEE',
                'JYOB': 'filiformis-JYOB',
                'JMJV': 'gaura-JMJV',
                'OEGL': 'glazioviana-OEGL',
                'MJHP': 'grandiflora-MJHP',
                'SEMK': 'grandis-SEMK',
                'OEHI': 'hirsuti-OEHI',
                'OEJA': 'jamesii-OEJA',
                'ICUS': 'laciniata-ICUS',
                'OELA': 'longissima-OELA',
                'HPNZ': 'longituba-HPNZ',
                'CJGZ': 'nana-CJGZ',
                'OENU': 'nutans-OENU',
                'OEOA': 'oakesiana-OEOA',
                'OEPA': 'parviflora-OEPA',
                'DZLN': 'picensis-DZLN',
                'YHLF': 'rhombipetala-YHLF',
                'FXJC': 'rosea-FXJC',
                'XZAQ': 'speciosa-XZAQ',
                'OEST': 'stuchii-OEST',
                'JKNQ': 'suffulta-JKNQ',
                'HKMQ': 'villaricae-HKMQ',
                'OEVS': 'villosa_strigosa-OEVS',
                'OEWO': 'wolfii-OEWO',
}

#print(species_map)

# path to trinotate functional annotation tables

trinotate_path = '/datahome/oenothera/transcripts/assemblies/snakemake-trinotate/results/reports'

	
report_list = [trinotate_path + '/' + key + '.trinotate_annotation_report.tsv' for key in species_map.keys()]


def is_comment(s):
    """ function to check if a line
         starts with some character.
         Here # for comment
    """
    # return true if a line starts with #
    return s.startswith('#')

# create function for finding phenolic proteins in the annotation reports
def search_report(report):

	ID = report.split('.')[0].split('/')[-1]
	taxon = species_map.get(ID)

	print("Finding annotation results for {}".format(taxon))

	phenolic_results = []

	with open(report) as anno_file:
		for line in dropwhile(is_comment, anno_file):
			#print(line)
			if line.split('\t')[2] != '.' and line.split('\t')[4] != '.':
				prot = line.split('\t')[4]
				sprot_result = str(line.split('\t')[2].split('Full=')[1].split(';')[0])
				result = ''.join(pattern.findall(sprot_result))
				if result:
					output = [taxon, prot, result, sprot_result]
					phenolic_results.append(output)
					#print(taxon, prot,result, sprot_result)
	return(phenolic_results)


# parallelize phenolics searches over all of the reports


with mp.Pool(processes=30) as pool:
	results = pool.starmap_async(search_report, zip(report_list))
	myresults = results.get()

pool.close()
pool.join() # Wait for all child processes to close.

#print(myresults)

phenolic_proteins = list(chain(*myresults))
print(phenolic_proteins)

df_phenolics = pd.DataFrame(phenolic_proteins, columns = ['Taxon', 'Protein ID', 'Search Term', 'Phenolic Protein'])

print(df_phenolics)

df_phenolics.to_csv('phenolic_proteins.csv', index = False)

#yivz_results = search_report(trinotate_path+ '/YIVZ.trinotate_annotation_report.tsv')
#print(yivz_results)
