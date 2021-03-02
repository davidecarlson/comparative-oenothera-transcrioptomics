#!/usr/bin/env python

import re
import sys
import glob
from itertools import chain
import pandas as pd
from collections import Counter
import csv
import multiprocessing as mp
from itertools import repeat

species_map = {
    'AJBK': 'elata-AJBK',
    'AQXA': 'laciniata-AQXA',
    'CHTD': 'grandis-CHTD',
    'CJGZ': 'nana-CJGZ',
    'DIBI': 'laciniata-DIBI',
    'DSVQ': 'biennis-DSVQ',
    'DWAP': 'biennis-DWAP',
    'DZLN': 'picensis-DZLN',
    'EBPG': 'grandiflora-EBPG',
    'ECSL': 'filiformis-ECSL',
    'EQYT': 'berlandieri-EQYT',
    'EWHG': 'grandis-EWHG',
    'EXGW': 'biennis-EXGW',
    'FDBS': 'biennis-FDBS',
    'FXJC': 'rosea-FXJC',
    'GBGW': 'speciosa-GBGW',
    'GVCB': 'affinis-GVCB',
    'GZUR': 'grandiflora-GZUR',
    'HKMQ': 'villaricae-HKMQ',
    'HOOU': 'grandiflora-HOOU',
    'HPNZ': 'longituba-HPNZ',
    'ICUS': 'laciniata-ICUS',
    'IDAU': 'elata_hookeri-IDAU',
    'IIFB': 'gaura-IIFB',
    'JKNQ': 'suffulta-JKNQ',
    'JMJV': 'gaura-JMJV',
    'JYOB': 'filiformis-JYOB',
    'KBRW': 'clelandii-KBRW',
    'KFAL': 'rosea-KFAL',
    'LACT': 'gaura-LACT',
    'LHDP': 'gaura-LHDP',
    'MJHP': 'grandiflora-MJHP',
    'MLUJ': 'biennis-MLUJ',
    'NQWS': 'filiformis-NQWS',
    'OEAR': 'argillicola-OEAR',
    'OEEE': 'elata-OEEE',
    'OEGL': 'glazioviana-OEGL',
    'OEHI': 'hirsuti-OEHI',
    'OEJA': 'jamesii-OEJA',
    'OELA': 'longissima-OELA',
    'OENU': 'nutans-OENU',
    'OEOA': 'oakesiana-OEOA',
    'OEPA': 'parviflora-OEPA',
    'OEST': 'stuchii-OEST',
    'OEVS': 'villosa_strigosa-OEVS',
    'OEVV': 'villosa_villosa-OEVV',
    'OEWO': 'wolfii-OEWO',
    'ONZG': 'laciniata-ONZG',
    'QBRF': 'laciniata-QBRF',
    'REZJ': 'grandiflora-REZJ',
    'RLZU': 'speciosa-RLZU',
    'ROLB': 'elata-ROLB',
    'SEMK': 'grandis-SEMK',
    'SJAN': 'serrulata-SJAN',
    'TAHF': 'grandis-TAHF',
    'TGOP': 'gaura-TGOP',
    'UASK': 'speciosa-UASK',
    'WPUV': 'grandis-WPUV',
    'XSNO': 'rosea-XSNO',
    'XZAQ': 'speciosa-XZAQ',
    'XZPP': 'speciosa-XZPP',
    'YHLF': 'rhombipetala-YHLF',
    'YIVZ': 'rosea-YIVZ'
}

# path to trinotate functional annotation tables

trinotate_path = '/datahome/oenothera/transcripts/assemblies/snakemake-trinotate/results/reports'

report_list = glob.glob(trinotate_path + '/*.trinotate_annotation_report.tsv')
#report_list = glob.glob('*.trinotate_annotation_report_test.tsv')
print(report_list)

report = 'AJBK.trinotate_annotation_report_test.tsv'

# count how many times each GO annotation appears and get the proportion of proteins with annotation


def count_annot(annot_list, taxon, report):
    annot_counter = Counter(annot_list)
    with open(report, 'r') as file:
        num_prots = 0
        for line in file:
            if line.split('\t')[4] != '.':
                # print(line.split('\t')[4])
                num_prots += 1
        print(num_prots)

    results = []
    for key, value in annot_counter.items():
        count = value
        prop = float(count) / num_prots
        annot = key
        combined = [annot, taxon, count, prop]
        results.append(combined)
        #print(OG, count, prop)
    return(results)

# create function that will count the GO annotations for each transcriptome annotation report


def get_annot_results(report):

    ID = report.split('.')[0].split('/')[-1]
    print(ID)
    taxon = species_map.get(ID)

    print("Finding annotation results for " + taxon)

    # create lists that will be populated as search through the report

    CC_list = []
    MF_list = []
    BP_list = []

    with open(report) as anno_file:
        for line in anno_file:
            if line.split('\t')[13] != '.':
                # print(line)
                blastp_annot = line.split('\t')[13]
                annot_list = re.split(r'(GO:[0-9]+\^)', blastp_annot)
                CC = [i.split('^')[1].replace('`', '')
                      for i in annot_list if i.startswith('cellular_component')]
                MF = [i.split('^')[1].replace('`', '')
                      for i in annot_list if i.startswith('molecular_function')]
                BP = [i.split('^')[1].replace('`', '')
                      for i in annot_list if i.startswith('biological_process')]

                CC_list.append(CC)
                MF_list.append(MF)
                BP_list.append(BP)

    # flatten the list of lists
    CC_list = sorted(list(chain(*CC_list)))
    # print(CC_list)
    MF_list = sorted(list(chain(*MF_list)))
    BP_list = sorted(list(chain(*BP_list)))

    #  get count and proportions of each GO annotation in the report
    CC_results = count_annot(CC_list, taxon, report)
    MF_results = count_annot(MF_list, taxon, report)
    BP_results = count_annot(BP_list, taxon, report)

    return(CC_results, MF_results, BP_results)


# parallelize finding all the annotations over each of the reports


with mp.Pool(processes=40) as pool:
    results = pool.starmap_async(get_annot_results, zip(report_list))
    myresults = results.get()

pool.close()
pool.join()  # Wait for all child processes to close.

# split results by the 3 GO categories and flatten the lists
myCC = list(chain(*[result[0] for result in myresults if result[0]]))
myMF = list(chain(*[result[1] for result in myresults if result[1]]))
myBP = list(chain(*[result[2] for result in myresults if result[2]]))


print(myCC)

df_CC = pd.DataFrame(myCC, columns=[
                     'GO Annotation', 'Taxon', 'Count of Proteins with Annotation', 'Proportion of Proteins with Annotation'])
df_CC.sort_values(
    by=['Taxon', 'Count of Proteins with Annotation'], ascending=False, inplace=True)
print(df_CC)

df_MF = pd.DataFrame(myMF, columns=[
                     'GO Annotation', 'Taxon', 'Count of Proteins with Annotation', 'Proportion of Proteins with Annotation'])
df_MF.sort_values(
    by=['Taxon', 'Count of Proteins with Annotation'], ascending=False, inplace=True)
print(df_MF)

df_BP = pd.DataFrame(myBP, columns=[
                     'GO Annotation', 'Taxon', 'Count of Proteins with Annotation', 'Proportion of Proteins with Annotation'])
df_BP.sort_values(
    by=['Taxon', 'Count of Proteins with Annotation'], ascending=False, inplace=True)
print(df_BP)

df_CC.to_csv('cellular_component_byTaxon.csv', index=False)
df_MF.to_csv('molecular_function_byTaxon.csv', index=False)
df_BP.to_csv('biological_process_byTaxon.csv', index=False)
