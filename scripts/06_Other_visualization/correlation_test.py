import matplotlib_venn as vplt
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import argparse
import sys
from statistics import mean
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns
import scipy.stats as stats
# setup
parser = argparse.ArgumentParser(
	description='This program will make a venn diagram of DEGs from single cell data')
# required arguments
parser.add_argument('--csv1', required=True, type=str,
	metavar='<str>', help='path to CSV containing mouse DEGs ending with .csv')
parser.add_argument('--csv2', required=True, type=str,
	metavar='<str>', help='path to CSV containing human DEGs ending with .csv')
parser.add_argument('--pdf', required=True, type=str,
	metavar='<str>', help='path to pdf output file ending with .pdf')
parser.add_argument('--pval', required=False, type=float, default=0.05,
	metavar='<float>', help='p-value cuttoff, default value is set to 0.05')

# finalization
arg = parser.parse_args()

##############################
## load complete data frame ##
##############################

complete_df1 = pd.read_csv(arg.csv1)
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/total_genes/Limma/all_cort_all_dataframe.csv
complete_df2 = pd.read_csv(arg.csv2)

#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/human_masterDEG.csv

#select significant genes and parse out males from females
complete_df1 = complete_df1.sort_values(by='pv', ascending=True)
sig_df1 = complete_df1[complete_df1['pv'] < arg.pval]
sig_df1 = sig_df1[sig_df1['sex'] == 'F']
sig_df1 = sig_df1[sig_df1['timepoint'] != 'E18']
print(sig_df1)

mouse_logfc = {}
for ind,row in sig_df1.iterrows():
	if row['gene'].upper() not in mouse_logfc:
		mouse_logfc[row['gene'].upper()] = []
	mouse_logfc[row['gene'].upper()].append(float(row['logfc']))

for gene,fc in mouse_logfc.items():
	mean = np.mean(np.array(fc))
	mouse_logfc[gene] = mean

#select significant genes and parse out males from females
complete_df2 = complete_df2.sort_values(by='adj.P.Val', ascending=True)
sig_df2 = complete_df2[complete_df2['adj.P.Val'] < arg.pval]
#sig_df2 = sig_df2[sig_df2['logFC'].abs() > arg.logFC]
human_logfc = {}
for ind,row in sig_df2.iterrows():
	if row['SYMBOL'] not in human_logfc:
		human_logfc[row['SYMBOL']] = []
	human_logfc[row['SYMBOL']].append(float(row['logFC']))

for gene,fc in human_logfc.items():
	mean = np.mean(np.array(fc))
	human_logfc[gene] = mean
	
a = []
b = []
for gene in human_logfc:
	if gene in mouse_logfc:
		a.append(human_logfc[gene])
		b.append(mouse_logfc[gene])
		
g = sns.jointplot(x=a, y=b, kind='reg', color='royalblue')
r, p = stats.spearmanr(a, b)
print(p)
g.ax_joint.annotate(f'$\\rho = {r:.3f}, p = {p:.4g}$',
                    xy=(0.1, 0.9), xycoords='axes fraction',
                    ha='left', va='center',
                    bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
g.ax_joint.scatter(a, b)
g.set_axis_labels(xlabel='human DEGs', ylabel='mouse DEGs', size=15)
plt.tight_layout()
plt.savefig(arg.pdf) 

	