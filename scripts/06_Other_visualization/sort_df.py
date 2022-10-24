import sys
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1])
df['Sum'] = df['Mouse_logFC'] + df['Human_logFC']
df1 = df.sort_values(['Sum', 'Human_logFC', 'Mouse_logFC'], ascending=[False, False, False])
top10 = df1[:20]
bot10 = df1[-20:]

top = pd.concat([top10, bot10])
print(top)

fig, axes = plt.subplots(ncols=2, sharey=True)
axes[0].barh(top['Human_Gene'].tolist(), top['Human_logFC'].tolist(), align='center', color=(top['Human_logFC']>0).map({True:'crimson', False: 'mediumblue'}))
axes[1].barh(top['Human_Gene'].tolist(), top['Mouse_logFC'].tolist(), align='center', color=(top['Mouse_logFC']>0).map({True:'crimson', False: 'mediumblue'}))
#axes[0].invert_xaxis()
axes[0].set(title='Top Human DEG logFC')
axes[1].set(title='Top Mouse DEG logFC')
axes[0].set(yticks=top['Human_Gene'].tolist())
axes[0].yaxis.tick_right()

for ax in axes.flat:
    ax.margins(0.03)
    ax.grid(True)

fig.tight_layout()
fig.subplots_adjust(wspace=0.35)
plt.savefig('top_DEG_human_mouse.pdf')

'''
#ran the following line
python3 /Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/06_Other_visualization/sort_df.py /Users/osman/Documents/GitHub/mecp2_chip_scrna_correlation/06_correlation_analysis/human_mouse.csv > test.csv
'''