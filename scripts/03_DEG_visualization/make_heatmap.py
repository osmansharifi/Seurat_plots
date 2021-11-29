#!/usr/bin/env python3

import argparse
import csv
import json
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
import mpl_toolkits
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import pandas as pd
import re
import sys

from natsort import natsorted
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

SMALL_SIZE = 9
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)
parser = argparse.ArgumentParser(
	description='DEG visualization with significance filtering',
	formatter_class=formatter)
parser.add_argument(
	'--frame', '-f', required=True, type=str, metavar='<path>', 
	help='path to dataframe for all the degs to plot')
parser.add_argument(
	'--celltypes', '-c', required=False, nargs='+', type=str,
	help='celltypes to plot')
parser.add_argument(
	'--sex', '-s', required=False, type=str, metavar='<str>',
	help='desired sex to filter by (M/F)')
parser.add_argument(
	'--timepoint', '-t', required=False, type=str, nargs='+',
	help='timepoints to filter by')
parser.add_argument(
	'--region', '-r', required=False, type=str, nargs='+',
	help='brain regions to filter by')
parser.add_argument(
	'--cutoff', '-u', required=False, type=float, metavar='<float>', default=0.05, 
	help='significance threshold for adjusted p-value')
parser.add_argument(
	'--num', '-n', required=False, type=int, default=50, metavar='<int>',
	help='number of significant DEGs to take from each celltypes')
parser.add_argument(
	'--fill', '-i', action='store_true',
	help='fill in plot with actual logFC value on places where gene not in celltype, otherwise fill with 0')
parser.add_argument(
	'--sort', '-o', required=False, type=str, default='max', metavar='<str>',
	help='sort method for significant genes, either max or std')
parser.add_argument(
	'--figdir', '-g', required=True, type=str, default=None, metavar='<path>',
	help='path to directory to save plot to')
parser.add_argument(
	'--save', '-v', required=False, type=str, default=None, metavar='<path>',
	help='save location for json file containing significant genes plotted, in a pandas dataframe')
parser.add_argument(
	'--pdf', required=True, type=str, default=None, metavar='<str>',
	help='name to save pdf plot')
parser.add_argument(
	'--verbose', '-b', required=False, type=int, default=0, metavar='<int>',
	help='verbosity level, 0 for no progress messages, 1 for progress messages')

def celltype_name(row, arg):
	ct = ''
	ct = row['celltype']
	if arg.sex == None: ct += '_'+row['sex']
	if arg.timepoint is not None:
		if len(arg.timepoint) > 1: ct += '_'+row['timepoint']
	elif arg.timepoint is None:
		ct += '_'+row['timepoint']
	if arg.region:
		if len(arg.region) > 1: ct += '_'+row['region']
	print(ct)
	return ct

arg = parser.parse_args()

df = pd.read_csv(arg.frame)
#df.drop(labels = 'Unnamed: 0')
print(df.head(5))
print(df.columns)
print(df.shape)

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

if ~newdf.empty:
	df = newdf
df.head(5)

df['ct_name'] = df.apply(lambda x: celltype_name(x, arg), axis=1)
print(df.tail(5))
print(df.columns)
print(df.shape)

sig_df = df.loc[df['pv'] <= arg.cutoff]
print(sig_df.head(20))
print(sig_df.shape)

sig_df['abslogfc'] = sig_df['logfc'].abs()

sig_df.sort_values(inplace=True,by=['celltype', 'timepoint', 'sex', 'region', 'abslogfc'], ascending=[True, True, True, True, False])
print()
print(sig_df.head(50))
print(sig_df.shape)

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

print(json.dumps(sig_genes, indent=2))
print(len(list(sig_genes.keys())))

celltypes = dict()
for sg in sig_genes:
	for ct in sig_genes[sg]:
		if ct not in celltypes: celltypes[ct] = True

print(json.dumps(celltypes, indent=2))

for ct in natsorted(celltypes.keys()):
	print(ct)

num_celltypes = len(list(celltypes.keys()))
num_genes = len(list(sig_genes.keys()))
df = df.set_index(['gene','ct_name'])
# print(df.head(5))
# print(df.columns)
# sys.exit()
print()
map   = np.zeros((num_genes, num_celltypes))
pvmap = np.zeros((num_genes, num_celltypes))
for i, siggene in enumerate(sig_genes.keys()):
	print(siggene)
	#sys.exit()
	for j, celltype in enumerate(natsorted(celltypes.keys())):
		print(siggene, celltype)
		try:
			down = df.loc[(siggene, celltype)]
			#print(down)
		except:
			down = pd.DataFrame()
			print('not down')
		print()
		if down.empty: pvmap[i, j] = 10.0
		else:          pvmap[i, j] = down.pv
		if celltype in sig_genes[siggene]:
			map[i, j] = sig_genes[siggene][celltype]
		else:
			if not arg.fill: map[i, j] = 0.0
			else:
				if down.empty: map[i, j] = 0.0
				else:          map[i, j] = down.logfc
print(map)
print(pvmap)
#sys.exit()
# cluster order the row labels because clustering moves rows
Z_rows = linkage(map, 'ward')

dn_rows = dendrogram(Z_rows)
	
gene_names = [nm for nm in sig_genes.keys()]
celltype_names = list(natsorted(celltypes.keys()))

gene_labels = [gene_names[ind] for ind in dn_rows['leaves']]

newpv_map = np.zeros((num_genes, num_celltypes))
frame = []
dic   = {}
for i, rg in enumerate(dn_rows['leaves']):
	dic = {}
	dic['gene'] = gene_names[rg]
	for j, ct in enumerate(celltype_names):
		dic[ct] = (map[rg, j], pvmap[rg, j])
		newpv_map[i, j] = pvmap[rg, j]
	frame.append(dic)

plot_df = pd.DataFrame(frame)
print(plot_df)

#sys.exit()
bool_map = np.zeros((num_genes, num_celltypes))
bool_map[np.where(newpv_map >= 0.05)] = 1

# plotting
fontsize_pt = plt.rcParams['figure.titlesize']
dpi = float(200)

matrix_height_pt = fontsize_pt * 2.5 * map.shape[0]
matrix_height_in = matrix_height_pt / dpi
matrix_width_pt = fontsize_pt * 3.0 * map.shape[1]
matrix_width_in = matrix_width_pt / dpi
print(matrix_height_in, matrix_width_in)
print(map.shape)
#sys.exit()
#top_margin = 0.02
#bottom_margin = 0.02

figure_height = matrix_height_in / (1 - 0.2)
figure_width = matrix_width_in / (1 - 0.5) 

plt.rc('font', size=matrix_height_in * dpi * 0.01)
plt.rc('figure', titlesize=matrix_width_in * dpi * 0.02)

plt.tick_params(
	axis='both',
	which='major',
	labelbottom=True,
	bottom=True,
	top=True,
	labeltop=True)

clm = sns.clustermap(
	map,
	method='ward',
	figsize=(figure_width,figure_height),
	row_cluster=True,
	col_cluster=False,
	row_linkage=Z_rows,
	z_score=None,
	standard_scale=None,
	norm=colors.TwoSlopeNorm(vmin=np.amin(map), vcenter=0, vmax=np.amax(map)),
	cmap=sns.color_palette('vlag', as_cmap=True),
	xticklabels=celltype_names,
	yticklabels=gene_labels,
	dendrogram_ratio=(0.10,0.0))

title_str = f'Top {arg.num} genes'

clm.fig.suptitle(title_str)
clm.fig.subplots_adjust(right=0.50, top=0.90)
plt.yticks(rotation=0)
clm.ax_heatmap.set_yticklabels(gene_labels)
clm.ax_heatmap.set_xticklabels(celltype_names)
clm.ax_cbar.set_position([0.90, 0.40, 0.50*clm.ax_row_dendrogram.get_position().width, 0.30])
clm.ax_cbar.set_title('logFC')
clm.ax_col_dendrogram.set_visible(False)

hm = clm.ax_heatmap

for i in range(bool_map.shape[0]):
	for j in range(bool_map.shape[1]):
		if bool_map[i, j] == 0.:
			hm.add_patch(Rectangle((j, i), 1, 1, edgecolor='black', fill=False, lw=1))
plt.close(1)
#plt.show()
plt.savefig(arg.pdf)
plt.close()

"""
1. to iterate through all cell types
	only grab top n genes from each
celltype_sex_time_region
"""















