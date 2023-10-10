import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# load datasets created from genes_in_neoloop.py
overlapping = pd.read_csv('neoloop_gene_overlapping.txt',sep='\t',header=None)
flanking = pd.read_csv('neoloop_gene_flanking.txt',sep='\t',header=None)

# overlapping
overlapping_unique = overlapping.apply(lambda x: x.drop_duplicates())
overlapping_unique = overlapping_unique.dropna(how='all')
overlapping_unique = overlapping_unique.drop(overlapping_unique.index[0])
overlapping_counts = overlapping_unique.stack().value_counts()
overlapping_counts_df = pd.DataFrame({'Value': overlapping_counts.index, 'Count': overlapping_counts.values})

# flanking
flanking_unique = flanking.apply(lambda x: x.drop_duplicates())
flanking_unique = flanking_unique.dropna(how='all')
flanking_unique = flanking_unique.drop(flanking_unique.index[0])
flanking_counts = flanking_unique.stack().value_counts()
flanking_counts_df = pd.DataFrame({'Value': flanking_counts.index, 'Count': flanking_counts.values})

overlapping_counts_df
merged_df = pd.merge(overlapping_counts_df, flanking_counts_df, on='Value', how='outer')
merged_df = merged_df.rename(columns={'Value': 'Gene', 'Count_x': 'Overlapping', 'Count_y': 'Flanking'})
merged_df['Total'] = merged_df['Overlapping'] + merged_df['Flanking']
merged_df = merged_df.fillna(0)
merged_df = merged_df.sort_values('Total', ascending = False)

merged_df.iloc[:20, :]

oncogenes_file = pd.read_csv('data/CancerGeneCensus.tsv',sep='\t',header=None)

########################################################################################
# test to make sure it works
# create a boolean mask by comparing the values in the desired column
mask = oncogenes_file[0] == 'CCR7'

# apply the mask to the DataFrame to get the filtered rows
filtered_rows = oncogenes_file[mask]
filtered_rows

########################################################################################
mydict = {} # create dict to make bar plot
mydict['Overlapping'] = merged_df.iloc[:20, 1].tolist()
mydict['Flanking'] = merged_df.iloc[:20, 2].tolist()

selectgenes = merged_df.iloc[:20, 0].tolist()

x = np.arange(len(selectgenes))  # the label locations
width = 0.3  # the width of the bars
multiplier = 0
fontsize = 16

colors = ['#94ACBF', '#4A6274']

fig, ax = plt.subplots(layout='constrained', figsize=(14, 8))

for i, (gene_type, genecounts) in enumerate(mydict.items()):
    offset = width * multiplier
    rects = ax.bar(x + offset, genecounts, width, label=gene_type, color=colors[i % len(colors)])  # Use colors[i % len(colors)] to cycle through the defined colors
    ax.bar_label(rects, padding=3, fontsize=fontsize-2)
    multiplier += 1

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks(x + width / 2)
ax.set_xticklabels(selectgenes, rotation=90, fontsize=fontsize)
ax.set_yticks([])
ax.tick_params(axis='y', labelsize=fontsize)

ax.set_title('Identification of Oncogenes in Neoloops across 79 GBM Samples', fontsize=fontsize+4)
ax.set_ylabel('Number of Samples', fontsize=fontsize+2)

ax.legend(loc='upper right', ncols=3, fontsize=fontsize)
ax.set_ylim(0, 24)

fig.set_facecolor('white')

plt.show()
