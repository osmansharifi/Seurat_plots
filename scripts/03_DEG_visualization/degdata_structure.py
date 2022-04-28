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

def read_dataexcel(excelpath, meta, dblist):
	print(excelpath)
	df = pd.read_excel(excelpath,header=0)
	df = df.drop(columns=['FC','Type','State','norm_foldChange'])
	print(df)
	print(df.shape)
	df[~df.SYMBOL.str.contains('|'.join(["^mt-"]))]
	print(df.shape)
	sys.exit()

arg = parser.parse_args()
assert(os.path.isdir(arg.dir))

data = []
dic = {}
for root, dirs, files in os.walk(arg.dir):
	if len(files) == 0: continue
	for file in files:
		if not file.endswith('.xlsx'): continue
		print(file)
		print(root)
		meta = root.split('/')
		tp_reg = meta[7]
		method = meta[8]
		ct = meta[9]
		
		sex = tp_reg.split('_')[4]
		tp  = tp_reg.split('_')[5]
		reg = tp_reg.split('_')[6]
		
		metadic = {
			'celltype': ct,
			'sex': sex,
			'timepoint': tp,
			'region': reg,
			'method': method
		}
		
		print(json.dumps(metadic,indent=2))
		
		frame = read_dataexcel(os.path.join(root,file), metadic, data)
		
		sys.exit()
	continue
	sys.exit()
	if not files.endswith('.xlsx'): continue
	print(file)
	continue
	assert(file.endswith('Limma_DEG.csv'))
	
	celltype = parse_celltype(file)
	
	splitted = file.split('_M_MUT_and_WT_')
	assert(len(splitted) == 2)
	
	#sex = splitted[0].split('_')[-1]
	rem_meta = splitted[1].split('_')
	sex = rem_meta[0]
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
sys.exit()
df = pd.DataFrame(data)
print(df.head(5))
print(df.columns)
print(df.shape)

df.to_csv(arg.output)

"""
all_data/
	all csvs
	celltype_M/F_MUT_and_WT_M/F_timepoint_region_Limma_DEG.csv
	sex, celltype, timepoint, region, gene, logfc, pv
	example of file name = Astro_M_MUT_and_WT_M_E18_WB_Limma_DEG.csv
	
directory_structure = /Users/osman/Desktop/LaSalle_lab/Rett_Data/Differential_expression/{Sex}/{timepoint and region}/{DEG_tool}/{celltype}/DEGs.xlsx
"""