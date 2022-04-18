#!/usr/bin/env python3

import argparse
import csv
import json
import os
import pandas as pd
import re
import sys

from natsort import natsorted
import numpy as np

formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=30)
parser = argparse.ArgumentParser(
	description='DEG visualization with significance filtering',
	formatter_class=formatter)
parser.add_argument(
	'--frame', '-f', required=True, type=str, metavar='<path>', 
	help='path to dataframe for all the degs to plot')
parser.add_argument(
	'--celltypes', '-c', required=False, nargs='+', type=str,
	default=False, help='celltypes to plot')
parser.add_argument(
	'--sex', '-s', required=False, type=str, metavar='<str>',
	default=False, help='desired sex to filter by (M/F)')
parser.add_argument(
	'--timepoint', '-t', required=False, type=str, nargs='+',
	default=False, help='timepoints to filter by')
parser.add_argument(
	'--region', '-r', required=False, type=str, nargs='+',
	default=False, help='brain regions to filter by')
parser.add_argument(
	'--cutoff', '-u', required=False, type=float, metavar='<float>',
	default=0.05, help='significance threshold for adjusted p-value')
parser.add_argument(
	'--num', '-n', required=False, type=int, default=50, metavar='<int>',
	help='number of significant DEGs to take from each celltypes')
parser.add_argument(
	'--fill', '-i', action='store_true',
	help=''.join((
		('fill in plot with actual logFC value on places where gene is not in'),
		('celltype, otherwise fill with 0'))))
parser.add_argument(
	'--save', '-v', required=True, type=str, metavar='<path>',
	help=''.join((
		('save location for csv files containing significant genes plotted.'),
		('csv is the serialized pandas dataframe after significance filter'))))
parser.add_argument(
	'--pdf', '-d', required=True, type=str, metavar='<str>',
	help='name and path to save pdf plot')
parser.add_argument('--verbose', '-b', required=False, type=int, default=0,
	metavar='<int>',
	help='verbosity level, 0 for no progress message, 1 for progress message')
parser.add_argument('--rotate', '-o', action='store_true',
	help='whether to rotate heatmap plot by -90 degrees in R')
parser.add_argument('--title', '-e', required=True, type=str,
	nargs='+', help='title for plot')

def celltype_name(row, arg):
	ct = ''
	ct = row['celltype']
	ct += '_'+row['sex']
	ct += '_'+row['timepoint']
	ct += '_'+row['region']
	return ct

arg = parser.parse_args()

save_path = arg.pdf.split('/')
save_dir = '/'.join(save_path[:-1])
save_name = save_path[-1].split('.')[0]
assert(os.path.isdir(save_dir))

if arg.save.endswith('.csv'):
	print('saving of frames need to go to a directory, not a file')
	sys.exit()

df = pd.read_csv(arg.frame)
df = df.drop(labels = 'Unnamed: 0', axis=1)
print("input dataframe")
print(df.head(5))
print(df.columns)
print(df.shape)
print()

newdf = pd.DataFrame()
if arg.celltypes:
	newdf = df.loc[df['celltype'].isin(arg.celltypes)]

if arg.sex:
	if newdf.empty: newdf = df.loc[df['sex'] == arg.sex]
	else:           newdf = newdf.loc[newdf['sex'] == arg.sex]

if arg.timepoint:
	if newdf.empty: newdf = df.loc[df['timepoint'].isin(arg.timepoint)]
	else:           newdf = newdf.loc[newdf['timepoint'].isin(arg.timepoint)]

if arg.region:
	if newdf.empty: newdf = df.loc[df['region'].isin(arg.region)]
	else:			newdf = newdf.loc[newdf['region'].isin(arg.region)]

if not newdf.empty:
	print('not empty')
	df = newdf

print("after filtering for metadata dataframe")
print(df.head(5))
print(df.columns)
print(df.shape)
print()

df['ct_name'] = df.apply(lambda x: celltype_name(x, arg), axis=1)
sig_df = df.loc[df['pv'] <= arg.cutoff].copy()
sig_df['abslogfc'] = sig_df['logfc'].abs()

sig_df.sort_values(
	inplace=True,
	by=[
		'celltype',
		'timepoint',
		'sex',
		'region',
		'abslogfc'],
	ascending=[True, True, True, True, False])

print("after pv significance filtering")
print(sig_df.head(5))
print(sig_df.shape)
print()

sig_genes = dict()
celltype = None
counter = 0
for indx, row in sig_df.iterrows():
	if celltype != row['ct_name']:
		celltype = row['ct_name']
		counter = 0
	if counter >= arg.num: 
		continue
	
	if row['gene'] not in sig_genes: sig_genes[row['gene']] = dict()
	sig_genes[row['gene']][row['ct_name']] = row['logfc']
	counter += 1

celltypes = dict()
for sg in sig_genes:
	for ct in sig_genes[sg]:
		if ct not in celltypes: celltypes[ct] = True

df = df.set_index(['gene','ct_name'])
filtered_logfc = []
filtered_pv    = []
for siggene in sig_genes:
	diclog = dict()
	dicpv  = dict()
	diclog['gene'] = siggene
	dicpv['gene']  = siggene
	for celltype in natsorted(celltypes.keys()):
		assert(celltype not in diclog)
		assert(celltype not in dicpv)
		diclog[celltype] = None
		dicpv[celltype]  = None
		
		try:
			down = df.loc[(siggene, celltype)]
			dicpv[celltype] = down.pv
			diclog[celltype] = down.logfc
		except:
			down = pd.DataFrame()
			dicpv[celltype] = 0.0
			diclog[celltype]   = 1.0
		
	filtered_logfc.append(diclog)
	filtered_pv.append(dicpv)

fdf_logfc = pd.DataFrame(filtered_logfc)
print(fdf_logfc.head(5))
print(fdf_logfc.columns)
print(fdf_logfc.shape)
fdf_logfc.to_csv('tmp_sn_logfc.csv')

fdf_pv = pd.DataFrame(filtered_pv)
print(fdf_pv.head(5))
print(fdf_pv.columns)
print(fdf_pv.shape)
fdf_pv.to_csv('tmp_sn_pv.csv')

cmd = ''.join((
	f"r_heatmaps.R -l tmp_sn_logfc.csv -p tmp_sn_pv.csv -s {arg.pdf}",
	f" -t {' '.join(arg.title)}"))

if arg.rotate:
	cmd += ' -r'

print(cmd)
os.system(cmd)
os.remove("tmp_sn_logfc.csv")
os.remove("tmp_sn_pv.csv")

fdf_logfc.to_csv(f"{arg.save}_logfc.csv")
fdf_pv.to_csv(f"{arg.save}_pv.csv")
