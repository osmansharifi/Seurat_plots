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

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=52)
parser = argparse.ArgumentParser(
	description='DEG visualization with significance filtering',
	formatter_class=formatter)
parser.add_argument(
	'--dir', '-d', required=True, type=str, metavar='<path>', 
	help='path to directory containing DEGs for celltypes within a condition')
parser.add_argument(
	'--cutoff', '-c', required=False, type=float, metavar='<float>', default=0.05, 
	help='significance threshold for adjusted p-value')
parser.add_argument(
	'--num', '-n', required=False, type=int, default=50, metavar='<int>',
	help='number of significant DEGs to plot across celltypes')
parser.add_argument(
	'--method', '-m', required=False, type=str, default='ward', metavar='<str>',
	help='clustering method')
parser.add_argument(
	'--fill', '-i', action='store_true',
	help='fill in plot with actual logFC value on places where gene not in celltype, otherwise fill with 0')
parser.add_argument(
	'--sort', '-o', required=False, type=str, default='max', metavar='<str>',
	help='sort method for significant genes, either max or std')
parser.add_argument(
	'--figdir', '-f', required=True, type=str, default=None, metavar='<path>',
	help='path to directory to save plots to')
parser.add_argument(
	'--save', '-s', required=False, type=str, default=None, metavar='<path>',
	help='save location for json file containing significant genes plotted, in a pandas dataframe')
parser.add_argument(
	'--verbose', '-v', required=False, type=int, default=0, metavar='<int>',
	help='verbosity level, 0 for no progress messages, 1 for progress messages')

arg = parser.parse_args()

assert(arg.sort == 'max' or arg.sort == 'std')
assert(os.path.isdir(arg.figdir))

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

def parse_celltype(csvname):
	parts = csvname.split('_')
	cell_type = ''
	for p in parts:
		if p == 'M' or p == 'F': return cell_type[:-1]
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

# plot per condition
for cond in data.keys():
	# significant genes
	significant_genes, total = significance_filtering(data[cond], arg.cutoff)
	
	if arg.verbose == 1:
		print(json.dumps(significant_genes,indent=2))
		print(len(list(total.keys())))
		print(len(list(significant_genes.keys())))
	
	# sort genes to find genes to plot
	if arg.sort == 'max':
		sigsorted = gene_maxsort(
			list(significant_genes.keys()),
			list(data[cond].keys()),
			significant_genes)
	else:
		sigsorted = gene_stdsort(
			list(significant_genes.keys()),
			list(data[cond].keys()),
			significant_genes)
	
	# sort cell types
	ctsorted = celltype_sort(
		list(significant_genes.keys()),
		list(data[cond].keys()),
		significant_genes)

	# create 2D matrix/map of data to eventually plot
	if arg.num == -1: num_genes = len(sigsorted)
	else:             num_genes = arg.num
	num_celltypes = len(ctsorted)
	
	map   = np.zeros((num_genes, num_celltypes))
	pvmap = np.zeros((num_genes, num_celltypes))
	for i, siggene in enumerate(sigsorted):
		if i == num_genes: break
		for j, celltype in enumerate(ctsorted):
			pvmap[i, j] = data[cond][celltype][siggene]['pv']
			if celltype in significant_genes[siggene]:
				map[i, j] = significant_genes[siggene][celltype]
			else:
				if not arg.fill: map[i, j] = 0.0
				else:
					if siggene in data[cond][celltype]:
						map[i, j] = data[cond][celltype][siggene]['logfc']
					else: map[i,j] = 0.0
	
	# cluster order the row/col labels because clustering moves rows+cols
	Z_rows = linkage(map, arg.method)
	Z_cols = linkage(np.transpose(map), arg.method)
	
	dn_rows = dendrogram(Z_rows)
	dn_cols = dendrogram(Z_cols)
	
	gene_names = [nm for nm in sigsorted.keys()]
	gene_names = gene_names[:num_genes]
	celltype_names = list(ctsorted.keys())
	
	gene_labels     = [gene_names[ind] for ind in dn_rows['leaves']]
	celltype_labels = [celltype_names[ind] for ind in dn_cols['leaves']]
	
	if arg.verbose == 1:
		print(dn_rows['leaves'])
		print(dn_cols['leaves'])
	
	newpv_map = np.zeros((num_genes, num_celltypes))
	frame = []
	dic   = {}
	for i, rg in enumerate(dn_rows['leaves']):
		dic = {}
		dic['gene'] = gene_names[rg]
		for j, ct in enumerate(dn_cols['leaves']):
			dic[celltype_names[ct]] = (map[rg, ct], pvmap[rg, ct])
			newpv_map[i, j] = pvmap[rg, ct]
		frame.append(dic)
	
	df = pd.DataFrame(frame)
	
	bool_map = np.zeros((num_genes, num_celltypes))
	bool_map[np.where(newpv_map >= 0.05)] = 1
	
	# plotting
	fontsize_pt = plt.rcParams['figure.titlesize']
	dpi = float(100)
	
	matrix_height_pt = fontsize_pt * 2.0 * map.shape[0]
	matrix_height_in = matrix_height_pt / dpi
	
	top_margin = 0.02
	bottom_margin = 0.02
	
	figure_height = matrix_height_in / (1 - top_margin - bottom_margin)
	
	plt.tick_params(
		axis='both',
		which='major',
		labelbottom=True,
		bottom=True,
		top=True,
		labeltop=True)
	
	clm = sns.clustermap(
		map,
		method=arg.method,
		figsize=(12,figure_height),
		row_cluster=True,
		col_cluster=True,
		row_linkage=Z_rows,
		col_linkage=Z_cols,
		z_score=None,
		standard_scale=None,
		norm=colors.TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5),
		cmap=sns.color_palette('vlag', as_cmap=True),
		xticklabels=celltype_labels,
		yticklabels=gene_labels,
		dendrogram_ratio=(0.20,0.0))
	
	title_str = f'Significant genes by cell type\n{cond}\nTop {arg.num} genes sorted by {arg.sort}'
	clm.fig.suptitle(title_str)
	clm.fig.subplots_adjust(right=0.80)
	plt.yticks(rotation=0)
	clm.ax_heatmap.set_yticklabels(gene_labels)
	clm.ax_heatmap.set_xticklabels(celltype_labels)
	clm.ax_cbar.set_position([0.90, 0.55, 0.25*clm.ax_row_dendrogram.get_position().width, 0.30])
	clm.ax_cbar.set_title('logFC')
	clm.ax_col_dendrogram.set_visible(False)
	
	hm = clm.ax_heatmap
	
	for i in range(bool_map.shape[0]):
		for j in range(bool_map.shape[1]):
			if bool_map[i, j] == 0.:
				hm.add_patch(Rectangle((j, i), 1, 1, edgecolor='black', fill=False, lw=1))
	
	
	plt.close(1)
	if arg.fill: fillstr = 'fill'
	else:        fillstr = ''
	filename = arg.figdir+cond+'_'+arg.sort+'_'+fillstr+'_'+str(arg.num)+'.pdf'
	plt.savefig(filename)
	plt.close()
	
	json_save = arg.save+'plotted_frame'+'_'+cond+'_'+arg.sort+'_'+fillstr+'_'+str(arg.num)+'.json'
	with open(json_save, 'w') as fw:
		df.to_json(fw)
