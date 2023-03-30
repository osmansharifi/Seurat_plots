import matplotlib_venn as vplt
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import argparse
import sys

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
parser.add_argument('--pval', required=False, type=float, default=0.1,
	metavar='<float>', help='p-value cuttoff, default value is set to 0.05')
parser.add_argument('--logFC', required=False, type=float, default=0,
	metavar='<float>', help='logFC cutoff, default value is set to 0.07')
# finalization
arg = parser.parse_args()

##############################
## load complete data frame ##
##############################

complete_df1 = pd.read_csv(arg.csv1)
#print(complete_df1)
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/total_genes/Limma/all_cort_all_dataframe.csv
complete_df2 = pd.read_csv(arg.csv2)
#print(complete_df2)
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/human_masterDEG.csv

#select significant genes and parse out males from females
complete_df1 = complete_df1.sort_values(by='pv', ascending=True)
sig_df1 = complete_df1[complete_df1['pv'] < arg.pval]
sig_df1 = sig_df1[sig_df1['logfc'].abs() > arg.logFC]
sig_df1 = sig_df1[sig_df1['timepoint'] != 'E18']
sig_df1 = sig_df1[sig_df1['sex'] == 'F']
print(len(sig_df1[sig_df1['timepoint'] == 'E18']))
print(len(sig_df1[sig_df1['timepoint'] == 'P30']))
print(len(sig_df1[sig_df1['timepoint'] == 'P60']))
print(len(sig_df1[sig_df1['timepoint'] == 'P150']))
print(len(sig_df1[sig_df1['timepoint'] == 'P120']))
print(len(sig_df1[sig_df1['sex'] == 'F']))
print(len(sig_df1[sig_df1['sex'] == 'M']))
#CSV2

#select significant genes and parse out males from females
complete_df2 = complete_df2.sort_values(by='adj.P.Val', ascending=True)
sig_df2 = complete_df2[complete_df2['adj.P.Val'] < arg.pval]
sig_df2 = sig_df2[sig_df2['logFC'].abs() > arg.logFC]
#print(len(sig_df))
mouse_genes = sig_df1[sig_df1['sex'] == 'F']
human_genes = sig_df2
print('There are', len(mouse_genes),'number of mouse DEGs.')
print('There are', len(human_genes),'number of human DEGs.')
#print(male_genes, female_genes)

print(len(mouse_genes), len(human_genes))
deduped_mouse = mouse_genes.drop_duplicates(subset=['gene'])
deduped_human = human_genes.drop_duplicates(subset=['SYMBOL'])
#print(len(deduped_male), len(deduped_female))
#print(deduped_male)
print('There are', len(deduped_mouse),'number of mouse DEGs after deduplication.')
print('There are', len(deduped_human),'number of human DEGs after deduplication.')

#plot a venn diagram of the differences between female and male genes
#plot a venn diagram of the differences between female and male genes
plt.figure(figsize=(12,12))
set1 = set(deduped_mouse['gene'])

newset1 = set()
for i in set1:
	newset1.add(i.upper()) 

set2 = set(deduped_human['SYMBOL'])
venn2([newset1, set2], ('mouse_DEGs', 'human_DEGs'), set_colors=('purple', 'skyblue'), alpha = 0.7)
plt.title("Venn diagram of postnatal DEGs from female mice cortices and human cortices")
plt.savefig(arg.pdf)  
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/venn_diagrams/all_female_male_venn.pdf'
#plt.show()

#Extract genes that are the intersection
print(len(set1.intersection(set2)))
print(set1.intersection(set2))

#check if a particular gene is present at the intersection
if 'Sst' in set1.intersection(set2):
    print('Gene found at the intersection')
else:
    print('Gene does not exist at the intersection')

#export a csv file containing the genes at the intersection
gene_list = list(newset1.intersection(set2))
with open('intersection_genes_postnatal.csv', 'w') as fp:
    fp.write('intersection genes for postnatal mouse and human' + '\n')
    for gene in gene_list:
        fp.write(gene + '\n')
