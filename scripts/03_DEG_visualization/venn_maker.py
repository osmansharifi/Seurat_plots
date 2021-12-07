import matplotlib_venn as vplt
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import argparse

# setup
parser = argparse.ArgumentParser(
	description='This will make a venn diagram of DEGs')
# required arguments
parser.add_argument('--csv', required=True, type=str,
	metavar='<str>', help='path to CSV containing DEGs ending with .csv')
parser.add_argument('--pdf', required=True, type=str,
	metavar='<str>', help='path to pdf output file ending with .pdf')
parser.add_argument('--pv', required=False, type=float, default=0.05,
	metavar='<float>', help='p-value cuttoff, default value is set to [%(default).2f]')
parser.add_argument('--logfc', required=False, type=float, default=0.7,
	metavar='<float>', help='logFC cutoff, default value is set to [%(default).1f]')
# finalization
arg = parser.parse_args()

#load complete data frame
complete_df = pd.read_csv(arg.csv)
print(complete_df)
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/total_genes/Limma/all_cort_all_dataframe.csv

#select significant genes and parse out males from females
complete_df = complete_df.sort_values(by='pv', ascending=True)
sig_df = complete_df[complete_df['pv'] < arg.pv]
sig_df = sig_df[sig_df['logfc'].abs() > arg.logfc]
#print(len(sig_df))
male_genes = sig_df[sig_df['sex'] == 'M']
female_genes = sig_df[sig_df['sex'] == 'F']
print('There are', len(male_genes),'number of male DEGs.')
print('There are', len(female_genes),'number of female DEGs.')
#print(male_genes, female_genes)

print(len(male_genes), len(female_genes))
deduped_male = male_genes.drop_duplicates(subset=['gene'])
deduped_female = female_genes.drop_duplicates(subset=['gene'])
#print(len(deduped_male), len(deduped_female))
#print(deduped_male)
print('There are', len(deduped_male),'number of male DEGs after deduplication.')
print('There are', len(deduped_female),'number of female DEGs after deduplication.')

#plot a venn diagram of the differences between female and male genes
#plot a venn diagram of the differences between female and male genes
plt.figure(figsize=(12,12))
set1 = set(deduped_male['gene'])
set2 = set(deduped_female['gene'])
venn2([set1, set2], ('male_DEGs', 'female_DEGs'), set_colors=('purple', 'skyblue'), alpha = 0.7)
plt.title("Venn diagram of all DEGs from males and females")
plt.savefig(arg.pdf)  
#'/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/venn_diagrams/all_female_male_venn.pdf'
#plt.show()

#Extract genes that are the intersection
print(len(set1.intersection(set2)))
print(set1.intersection(set2))

#check if a particular gene is present at the intersection
if 'Mbp' in set1.intersection(set2):
    print('Gene found at the intersection')
else:
    print('Gene does not exist at the intersection')

#export a csv file containing the genes at the intersection
gene_list = list(set1.intersection(set2))
with open('intersection_genes_test.csv', 'w') as fp:
    fp.write('intersection genes for males and females' + '\n')
    for gene in gene_list:
        fp.write(gene + '\n')
