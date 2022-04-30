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

def read_dataexcel(excelpath):
	df = pd.read_excel(excelpath,header=0)
	new = df[['SYMBOL', 'logFC', 'P.Value', 'adj.P.Val']].copy()
	del df
	for ind, row in new.iterrows():
		if re.match(r'^mt-.', row['SYMBOL']):
			new = new.drop(ind)
	
	return new

arg = parser.parse_args()
assert(os.path.isdir(arg.dir))

data = []
dic = {}
column_names = ['SYMBOL', 'logFC', 'P.Value', 'adj.P.Val', 'celltype', 'sex', 'timepoint', 'region', 'method']
total_df = pd.DataFrame(columns=column_names)

for root, dirs, files in os.walk(arg.dir):
	if len(files) == 0: continue
	for file in files:
		if not file.endswith('.xlsx'): continue
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
		
		frame = read_dataexcel(os.path.join(root,file))
		frame_rows = frame.shape[0]
		for k,v in metadic.items():
			new_list = [v] * frame_rows
			frame[k] = new_list
		
		total_df = pd.concat([total_df, frame])

total_df.to_csv(arg.output)
print(total_df)
print(total_df.shape)
print(total_df.columns)
