#!/usr/bin/env python3

import argparse
import csv
import json
import os
import pandas as pd
import re
import sys

formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)
parser = argparse.ArgumentParser(
	description='Make data structure for DEG data across timepoints and celltypes, sex, brain region',
	formatter_class=formatter)
parser.add_argument(
	'--dir', '-d', required=True, type=str, metavar='<path>', 
	help='path to directory containing DEG experiments')
parser.add_argument(
	'--output', '-o', required=True, type=str, metavar='<path>',
	help='path to location to save resulting dataframe as csv')

def parse_celltype(csvname):
	parts = csvname.split('_')
	cell_type = ''
	for p in parts:
		if p == 'M' or p == 'F': return cell_type[:-1]
		else:                    cell_type += p+'_'

def read_datacsv(csvpath, meta, dblist):
	dic = {}
	with open(csvpath, 'rt') as csvfile:
		csv_reader = csv.DictReader(csvfile)
		for row in csv_reader:
			if re.match(r'^mt-.', row['']): continue
			dic = {}
			dic = meta.copy()
			dic['gene'] = row['']
			dic['pv'] = row['adj.P.Val']
			dic['logfc'] = row['logFC']
			dblist.append(dic)
	return dblist

arg = parser.parse_args()
assert(os.path.isdir(arg.dir))

data = []
dic = {}
for file in os.listdir(arg.dir):
	dic = {}
	print(file)
	assert(file.endswith('Limma_DEG.csv'))
	
	celltype = parse_celltype(file)
	
	splitted = file.split('_MUT_and_WT_')
	assert(len(splitted) == 2)
	
	sex = splitted[0].split('_')[-1]
	rem_meta = splitted[1].split('_')
	timepoint = rem_meta[1]
	region = rem_meta[2]
	
	dic = {
		'celltype':celltype,
		'sex':sex,
		'timepoint':timepoint,
		'region':region
	}
	
	data = read_datacsv(os.path.join(arg.dir,file), dic, data)
	#print(json.dumps(data, indent=2))

df = pd.DataFrame(data)
print(df.head(5))
print(df.columns)
print(df.shape)

df.to_csv(arg.output)

# data = dict()
# for dir, subdirs, files in os.walk(path):
# 	if dir[-1] != '/': dir += '/'
# 	csvs = [f for f in files if 'csv' in f]
# 	if len(csvs) == 0: continue
# 	
# 	dirnames = dir.split('/')
# 	cond = dirnames[-2]
# 	if arg.verbose == 1:
# 		print(dirnames)
# 		print(cond)
# 	assert(len(subdirs) == 0)
# 	
# 	data[cond] = dict()
# 	for celltype in files:
# 		subpath = dir + '/' + celltype
# 		celltype_name = parse_celltype(celltype)
# 		data[cond][celltype_name] = dict()
# 		data[cond][celltype_name] = read_datacsv(subpath)
# 

"""
all_data/
	all csvs
	celltype_M/F_MUT_and_WT_M/F_timepoint_region_Limma_DEG.csv
	sex, celltype, timepoint, region, gene, logfc, pv
"""