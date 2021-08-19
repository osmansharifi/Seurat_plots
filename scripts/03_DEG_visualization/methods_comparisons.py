#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# Variables need to be changed to reflect data
# Names of the CSV files being input 
tab_names = ["M_MUT_and_WT_M_E18_WB_DEG_method_comparison",
            "M_MUT_and_WT_M_P30_CORT_DEG_method_comparison",
            "M_MUT_and_WT_M_P60_CORT_DEG_method_comparison",
            "M_MUT_and_WT_M_P120_CORT_DEG_method_comparison"]
# Metadata information corresponding to the CSV files (must be in the same order as tab_names)             
meta_data = ["Male, E18, Whole Brain",
            "Male, P30, Cortex",
            "Male, P60, Cortex",
            "Male, P120, Cortex"]

for tab, metadata in zip(tab_names, meta_data):
    df = pd.read_csv(f'../../DEG_data/{tab}.csv')
    # Make lists using Excel column values
    cell_types = df['Cluster'].tolist()
    deseq_total = df['DESeq2 (tot)'].tolist()
    limma_total = df['Limma (tot)'].tolist()
    edger_total = df['EdgeR (tot)'].tolist()
    deseq_only = df['DESeq2 Only'].tolist()
    limma_only = df['Limma Only'].tolist()
    edger_only = df['EdgeR Only'].tolist()
    deseq_and_limma = df['DESeq2 & Limma'].tolist()
    deseq_and_edger = df['DESeq2 & EdgeR'].tolist()
    limma_and_edger = df['Limma & EdgeR'].tolist()
    # Assign scatter plot lines
    one,=plt.plot(cell_types, deseq_total, label = 'DESeq2')
    two,=plt.plot(cell_types, limma_total, label = 'Limma')
    three,=plt.plot(cell_types, edger_total, label = 'EdgeR')
    # Make scatter plot for all methods
    fontP = FontProperties()
    fontP.set_size('medium')
    plt.title(f'Total DEGs Identified Per Method\n {metadata}')
    plt.xlabel('Cell Type (Cluster)')
    plt.xticks(rotation = 75)
    plt.ylabel('Number of DEGs')
    plt.subplots_adjust(bottom=0.275)
    plt.legend([one, two, three],
           ['DESeq2',
           'Limma',
           'EdgeR',],
           bbox_to_anchor=(1.05, 1),
           loc='upper left',
           prop=fontP)  
    plt.savefig(f'../../figures/methods_comparison/{tab}_total.pdf', bbox_inches = "tight", format = "pdf")
    plt.close()
    # Make scatter plot for unique DEGs (commands need to be rerun since plt.close() is used)
    df = pd.read_csv(f'../../DEG_data/{tab}.csv')
    # Make lists using Excel column values
    cell_types = df['Cluster'].tolist()
    deseq_total = df['DESeq2 (tot)'].tolist()
    limma_total = df['Limma (tot)'].tolist()
    edger_total = df['EdgeR (tot)'].tolist()
    deseq_only = df['DESeq2 Only'].tolist()
    limma_only = df['Limma Only'].tolist()
    edger_only = df['EdgeR Only'].tolist()
    deseq_and_limma = df['DESeq2 & Limma'].tolist()
    deseq_and_edger = df['DESeq2 & EdgeR'].tolist()
    limma_and_edger = df['Limma & EdgeR'].tolist()
    all_methods_total = df['All Methods'].tolist()
    # Assign scatter plot lines 
    four,=plt.plot(cell_types, deseq_only, label = 'DEGs Unique to DESeq2')
    five,=plt.plot(cell_types, limma_only, label = 'DEGs Unique to Limma')
    six,=plt.plot(cell_types, edger_only, label = 'DEGs Unique to EdgeR')
    seven,=plt.plot(cell_types, deseq_and_limma, label = 'DEGs Unique to DESeq2 and Limma')
    eight,=plt.plot(cell_types, deseq_and_edger, label = 'DEGs Unique to DESeq2 and EdgeR')
    nine,=plt.plot(cell_types, limma_and_edger, label = 'DEGs Unique to Limma and EdgeR')
    ten,=plt.plot(cell_types, all_methods_total, label = 'All Methods')
    # Make plot
    fontP = FontProperties()
    fontP.set_size('medium')
    plt.title(f'Comparison of Methods Used for DEG Identification\n {metadata}')
    plt.xlabel('Cell Type (Cluster)')
    plt.xticks(rotation = 75)
    plt.ylabel('Number of DEGs')
    plt.subplots_adjust(bottom=0.275)
    plt.legend([four, five, six, seven, eight, nine, ten],
           ['DEGs Unique to DESeq2', 
           'DEGs Unique to Limma',
           'DEGs Unique to EdgeR',
           'DEGs Unique to DESeq2 and Limma',
           'DEGs Unique to DESeq2 and EdgeR',
           'DEGs Unique to Limma and EdgeR',
           'All Methods'],
           bbox_to_anchor=(1.05, 1),
           loc='upper left',
           prop=fontP)  
    plt.savefig(f'../../figures/methods_comparison/{tab}_unique.pdf', bbox_inches = "tight", format = "pdf")
    plt.close()