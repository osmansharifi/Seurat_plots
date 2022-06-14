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
	'--degmethod', '-m', required=False, type=str, nargs='+',
	default=False, help='which deg method to filter by')
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
		('fill in plot with actual logFC value on places where SYMBOL is not in'),
		('celltype, otherwise fill with 0'))))
parser.add_argument(
	'--save', '-v', required=True, type=str, metavar='<path>',
	help=''.join((
		('save location for csv files containing significant SYMBOLs plotted.'),
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
	if row['method'] == 'limmaVoomCC': ct += '_'+'1'
	else:                              ct += '_'+'2'
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
df = df.drop(labels = 'Unnamed: 0', axis=1).dropna()
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

if arg.degmethod:
	if newdf.empty: newdf = df.loc[df['method'].isin(arg.degmethod)]
	else:           newdf = newdf.loc[newdf['method'].isin(arg.degmethod)]

if not newdf.empty:
	print('not empty')
	df = newdf

print("after filtering for metadata dataframe")
print(df.head(5))
print(df.columns)
print(df.shape)
print()

df['ct_name'] = df.apply(lambda x: celltype_name(x, arg), axis=1)
sig_df = df.loc[df['adj.P.Val'] <= arg.cutoff].copy()
sig_df['abslogFC'] = sig_df['logFC'].abs()

sig_df.sort_values(
	inplace=True,
	by=[
		'celltype',
		'timepoint',
		'sex',
		'region',
		'method',
		'abslogFC'],
	ascending=[True, True, True, True, True, False])

print("after pv significance filtering")
print(sig_df.head(5))
print(sig_df.shape)
print()

sig_SYMBOLs = dict()
celltype = None
counter = 0
for indx, row in sig_df.iterrows():
	if celltype != row['ct_name']:
		celltype = row['ct_name']
		counter = 0
	if counter >= arg.num: 
		continue
	
	if row['SYMBOL'] not in sig_SYMBOLs: sig_SYMBOLs[row['SYMBOL']] = dict()
	sig_SYMBOLs[row['SYMBOL']][row['ct_name']] = row['logFC']
	counter += 1

celltypes = dict()
for sg in sig_SYMBOLs:
	for ct in sig_SYMBOLs[sg]:
		if ct not in celltypes: celltypes[ct] = True

df = df.set_index(['SYMBOL','ct_name'])
print(df)
#sys.exit()
filtered_logFC = []
filtered_pv    = []
for sigSYMBOL in sig_SYMBOLs:
	diclog = dict()
	dicpv  = dict()
	diclog['SYMBOL'] = sigSYMBOL
	dicpv['SYMBOL']  = sigSYMBOL
	for celltype in natsorted(celltypes.keys()):
		assert(celltype not in diclog)
		assert(celltype not in dicpv)
		diclog[celltype] = None
		dicpv[celltype]  = None
		#print('here', sigSYMBOL, celltype)
		#down = df.loc[(sigSYMBOL, celltype)]
		#print(down)
		#dicpv[celltype] = down.
		try:
			down = df.loc[(sigSYMBOL, celltype)]
			dicpv[celltype]  = down['adj.P.Val']
			diclog[celltype] = down.logFC
			#print(down.logFC)
		except:
			down = pd.DataFrame()
			dicpv[celltype]  = 1.0
			diclog[celltype] = 0.0
		
	filtered_logFC.append(diclog)
	filtered_pv.append(dicpv)
#sys.exit()
fdf_logFC = pd.DataFrame(filtered_logFC)
print(fdf_logFC.head(25))
print(fdf_logFC.columns)
print(fdf_logFC.shape)
#sys.exit()
fdf_logFC.to_csv('tmp_sn_logFC.csv')

fdf_pv = pd.DataFrame(filtered_pv)
print(fdf_pv.head(5))
print(fdf_pv.columns)
print(fdf_pv.shape)
fdf_pv.to_csv('tmp_sn_pv.csv')

cmd = ''.join((
	f"/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/03_DEG_visualization/r_heatmaps.R -l tmp_sn_logFC.csv -p tmp_sn_pv.csv -s {arg.pdf}",
	f" -t {' '.join(arg.title)}"))

if arg.rotate:
	cmd += ' -r'

print(cmd)
os.system(cmd)
os.remove("tmp_sn_logFC.csv")
os.remove("tmp_sn_pv.csv")

fdf_logFC.to_csv(f"{arg.save}_logFC.csv")
fdf_pv.to_csv(f"{arg.save}_pv.csv")
