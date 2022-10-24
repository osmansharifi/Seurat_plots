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

df_human = pd.read_csv(arg.csv1)

#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/total_genes/Limma/all_cort_all_dataframe.csv
df_mouse = pd.read_csv(arg.csv2)
df_mouse.to_csv('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/mouse_GOterms.csv')
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/07_human_cell_labeling/human_rett_cort_filt/human_masterDEG.csv

#select significant genes and parse out males from females

df_mouse = df_mouse[df_mouse['Sex'] != 'males']
df_mouse = df_mouse[df_mouse['Time Point'] != 'E18']

#Remove the GO accession number from Terms
new_goterms = []
for go in df_human['Term']:
	new_goterms.append(' '.join(go.split(' ')[:-1]))
df_human['Term'] = new_goterms
print(len(df_mouse['Term']))
print(len(df_human['Term']))
print(df_human['Term'])
df_human.to_csv('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/human_GOTerms.csv')
plt.figure(figsize=(12,12))
set1 = set(df_mouse['Term'])
set2 = set(df_human['Term'])
venn2([set1, set2], ('mouse_GO_Terms', 'human_GO_Terms'), set_colors=('purple', 'skyblue'), alpha = 0.9)
plt.title("Venn diagram of postnatal GO Terms from female mice cortices and human cortices")
plt.savefig(arg.pdf)  
print(len(set1))
print(len(set2))

#Extract genes that are the intersection
print(len(set1.intersection(set2)))
print(set1.intersection(set2))
sys.exit()
#check if a particular gene is present at the intersection
if 'Sst' in set1.intersection(set2):
    print('Gene found at the intersection')
else:
    print('Gene does not exist at the intersection')

#export a csv file containing the genes at the intersection
gene_list = list(set1.intersection(set2))
with open('intersection_genes_postnatal.csv', 'w') as fp:
    fp.write('intersection genes for postnatal mouse and human' + '\n')
    for gene in gene_list:
        fp.write(gene + '\n')
