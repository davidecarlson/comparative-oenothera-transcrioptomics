#!/usr/bin/env python

import glob
import csv
import multiprocessing
import pandas as pd
import statistics as stats


def sum_quant(file):
	name = file.split('/')[9]
	df = pd.read_table(file, header=0, sep='\t')
	num_genes = len(df.index)
	expressed = df[df.TPM > 0]
	num_expressed = len(expressed.index)
	avg_tpm = stats.mean(df['TPM'])
	med_tpm = stats.median(df['TPM'])
	results = [name, num_genes, num_expressed, avg_tpm, med_tpm]

	return(results)

p = multiprocessing.Pool(25)

results = p.map_async(sum_quant, (file for file in glob.glob('/bowman/datahome_migrated/oenothera/transcripts/assemblies/snakemake-trinity/results/quantification/*/quant.sf.genes')))

all_results = results.get()

# could also use p.map and do the analysis synchronously and skip the results.get() step. Doesn't seem to make a time difference either way for this dataset

sorted_results = sorted(all_results, key=lambda x: x[0])

header = ['Name', 'Number of Genes', 'Number of Expressed Genes', 'Mean TPM', 'Median TPM']

with open("TPM_summary.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for sample in sorted_results:
    	writer.writerow(sample)

