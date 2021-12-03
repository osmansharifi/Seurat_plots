%matplotlib inline
import matplotlib_venn as vplt
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import pandas as pd

#load complete data frame
complete_df = pd.read_csv('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/total_genes/Limma/all_cort_all_dataframe.csv')
print(complete_df)

#select significant genes and parse out males from females
complete_df = complete_df.sort_values(by='pv', ascending=True)
sig_df = complete_df[complete_df['pv'] < 0.05]
sig_df = sig_df[sig_df['logfc'].abs() > 0.7]
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
plt.savefig('/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/venn_diagrams/all_female_male_venn.pdf')  
#plt.show()

#Extract genes that are the intersection
print(len(set1.intersection(set2)))
print(set1.intersection(set2))
