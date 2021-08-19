#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

tab_names = ["M_MUT_and_WT_M_E18_WB_DEG_method_comparison",
             "M_MUT_and_WT_M_P30_CORT_DEG_method_comparison",
             "M_MUT_and_WT_M_P60_CORT_DEG_method_comparison",
             "M_MUT_and_WT_M_P120_CORT_DEG_method_comparison"]

for tab in tab_names:
   df = pd.read_csv(f'../../DEG_data/{tab}.csv')
   print(tab)
   print(df)

'''
# Make lists using Excel column values
freq_A = df[' f(A)'].tolist()
A_rec_one = df['Delta f(A) Rec 0.1'].tolist()
A_dom_one = df['Delta f(A) Dom 0.1'].tolist()
A_add_one = df['Delta f(A) Add 0.1'].tolist()
A_rec_twofive = df['Delta f(A) Rec 0.25'].tolist()
A_dom_twofive = df['Delta f(A) Dom 0.25'].tolist()
A_add_twofive = df['Delta f(A) Add 0.25'].tolist()

# Delineate scatter plot lines
one,= plt.plot(freq_A, A_rec_one, label = 'Recessive Allele S = 0.1')
two,= plt.plot(freq_A, A_dom_one, label = 'Dominant Allele S = 0.1')
three,= plt.plot(freq_A, A_add_one, label = 'Additive Allele S = 0.1')
four,= plt.plot(freq_A, A_rec_twofive, label = 'Recessive Allele S = 0.25')
five,= plt.plot(freq_A, A_dom_twofive, label = 'Dominant Allele S = 0.25')
six,= plt.plot(freq_A, A_add_twofive, label = 'Additive Allele S = 0.25')

# Make scatter plot
fontP = FontProperties()
fontP.set_size('medium')
plt.title('Properties of Selection in Large Populations')
plt.xlabel('Frequency of the Advantageous Allele')
plt.ylabel('Change in Frequency of the Advantageous Allele')
plt.legend([one, two, three, four, five, six],
           ['Recessive Allele S = 0.1', 'Dominant Allele S = 0.1', 'Additive Allele S = 0.1', 'Recessive Allele S = 0.25', 'Dominant Allele S = 0.25', 'Additive Allele S = 0.25'],
           bbox_to_anchor=(1.05, 1),
           loc='upper left',
           prop=fontP)  
plt.show()
'''