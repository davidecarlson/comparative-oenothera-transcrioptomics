#!/usr/bin/env python

# get number of transcripts for each assembly that are longer than some cutoff

import argparse
import glob
import numpy as np
from Bio import SeqIO

parser = argparse.ArgumentParser(description='get number of transcripts for each assembly that are longer than some cutoff')

parser.add_argument('--dir', required=True, dest='dir', action='store')
parser.add_argument('--min_len', required=True, dest='min_len', action='store')

args=parser.parse_args()

asm_dir = args.dir
min_len = args.min_len

asm_list = glob.glob(asm_dir + '/*_trinity.Trinity.fasta')



def trans_count(asm):
	count = 0
	for record in SeqIO.parse(asm, "fasta"):
		trans_len = len(record.seq)
		if trans_len >= int(min_len):
			count+=1
	return(count)

#test = trans_count('/bowman/datahome_migrated/oenothera/transcripts/assemblies/snakemake-trinity/results/assemblies/AJBK_trinity.Trinity.fasta')

results = []

for asm in asm_list:
	asm_count = trans_count(asm)
	results.append(asm_count)

#print(results)

avg = np.round(np.mean(results), 2)
print(f'The average number of transcripts per assembly at least {min_len} in length is {avg}')
	
	
		
		
			



