#!/usr/bin/env python3

import argparse
import csv
import json
import os
import re
import sys

import numpy as np
import pandas as pd

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

parser = argparse.ArgumentParser(
	description='DEG visualization with significance filtering')
parser.add_argument(
	'--dir', '-d', required=True, type=str, metavar='<path>', 
	help='path to directory containing DEGs for celltypes within a condition')
parser.add_argument(
	'--cutoff', '-c', required=False, type=float, metavar='<float>', default=0.05, 
	help='significance threshold for adjusted p-value')
parser.add_argument(
	'--num', '-n', required=False, type=int, default=20, metavar='<int>',
	help='number of significant DEGs to plot across celltypes')
parser.add_argument(
	'--verbose', '-v', required=False, type=int, default=0, metavar='<int>',
	help='verbosity level, 0 for no progress messages, 1 for progress messages')

arg = parser.parse_args()

def parse_celltype(csvname):
	parts = csvname.split('_')
	cell_type = ''
	for p in parts:
		if p == 'M': return cell_type[:-1]
		else:        cell_type += p+'_'

def read_datacsv(csvpath):
	sigdic = dict()
	with open(csvpath, 'rt') as csvfile:
		csv_reader = csv.DictReader(csvfile)
		for row in csv_reader:
			if re.match(r'^mt-.', row['']): continue
			sigdic[row['']] = {
				'pv': float(row['adj.P.Val']),
				'logfc': float(row['logFC'])}
	return sigdic

def significance_filtering(data, threshold):
	siggenes_ct = dict()
	total = dict()
	
	for celltype in data:
		for gene in data[celltype]:
			if gene not in total: total[gene] = 0
			total[gene] += 1
			if data[celltype][gene]['pv'] < threshold:
				if gene not in siggenes_ct: siggenes_ct[gene] = dict()
				siggenes_ct[gene][celltype] = data[celltype][gene]['logfc']
	
	return siggenes_ct, total    

def gene_maxsort(genes, celltypes, data):
	genesorted = dict()
	for gene in genes:
		if gene not in data: continue
		logfcs = []
		for celltype in celltypes:
			if celltype not in data[gene]: logfcs.append(0.0)
			else:                          logfcs.append(np.abs(data[gene][celltype]))
		#print(logfcs)
		max = np.amax(np.array(logfcs))
		#print(var, gene)
		genesorted[gene] = max
	keys = {k:v for k,v in sorted(genesorted.items(), key=lambda x: x[1], reverse=True)}
	return keys

def gene_stdsort(genes, celltypes, data):
	genesorted = dict()
	for gene in genes:
		if gene not in data: continue
		logfcs = []
		for celltype in celltypes:
			if celltype not in data[gene]: logfcs.append(0.0)
			else:                          logfcs.append(data[gene][celltype])
		
		std = np.std(np.array(logfcs))
		genesorted[gene] = std
	keys = {k:v for k,v in sorted(genesorted.items(), key=lambda x: x[1], reverse=True)}
	return keys

def celltype_sort(genes, celltypes, data):
	ct_logfcs = dict()
	for celltype in celltypes:
		logfcs = []
		for gene in genes:
			if celltype not in data[gene]: logfcs.append(0.0)
			else:                          logfcs.append(np.abs(data[gene][celltype]))
		if len(logfcs) == 0: continue
		avg_logfc = np.mean(np.array(logfcs))
		if avg_logfc == 0.0: continue
		ct_logfcs[celltype] = avg_logfc
	
	keys = {k:v for k,v in sorted(ct_logfcs.items(), key=lambda x: x[1], reverse=True)}
	return keys

def make_frame(genes, celltypes, data):
	frame = list()
	for g in genes:
		dic  = {'gene':g}
		gset = []
		for ct in celltypes:
			if ct not in dic: dic[ct] = 0.0
			if g in data[ct]:
				dic[ct] = data[ct][g]['logfc']
				gset.append(data[ct][g]['logfc'])
			else: gset.append(0.0)
		
		dic['max'] = np.amax(np.abs(np.array(gset)))
		dic['std'] = np.std(np.array(gset))
		dic['avg'] = np.mean(np.array(gset))
		
		frame.append(dic)	
	
	df = pd.DataFrame(frame)
	return df

path = arg.dir

if path[-1] != '/': path += '/'
# collect all the data
data = dict()
for dir, subdirs, files in os.walk(path):
	if dir[-1] != '/': dir += '/'
	csvs = [f for f in files if 'csv' in f]
	if len(csvs) == 0: continue
	
	dirnames = dir.split('/')
	cond = dirnames[-2]
	if arg.verbose == 1:
		print(dirnames)
		print(cond)
	assert(len(subdirs) == 0)
	
	data[cond] = dict()
	for celltype in files:
		subpath = dir + '/' + celltype
		celltype_name = parse_celltype(celltype)
		data[cond][celltype_name] = dict()
		data[cond][celltype_name] = read_datacsv(subpath)

assert(len(list(data.keys())) == 1)

cond = list(data.keys())[0]

significant_genes, total = significance_filtering(data[cond], arg.cutoff)

ctsorted = celltype_sort(
	list(significant_genes.keys()),
	list(data[cond].keys()),
	significant_genes)

sig_df = make_frame(list(significant_genes.keys()), list(ctsorted.keys()), data[cond])
blankindex = [''] * len(sig_df)
sig_df.index = blankindex
print(sig_df.sort_values(by='std', ascending=True))
print(sig_df.columns)
