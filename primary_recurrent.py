import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.stats as stats
import os

# import and clean metadata
metadata = pd.read_excel('CATALYST - MS96 Selected.v5-ATAC.Seq.xlsx')
metadata = metadata[['HiC Seq Corrected', 'WGS Seq ID', 'RNA Seq ID', 'Tissue Type']]
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('Int64')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('str')
metadata['WGS Seq ID'] = pd.to_numeric(metadata['WGS Seq ID'], errors='coerce').astype('Int64')
metadata['WGS Seq ID'] = pd.to_numeric(metadata['WGS Seq ID'], errors='coerce').astype('str')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('Int64')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('str')

# rename metadata columns to merge later
metadata.rename(columns={'HiC Seq Corrected': 'HiC_ID', 
                         'WGS Seq ID': 'WGS_ID',
                        'RNA Seq ID': 'RNA_ID'}, inplace=True)

# import and clean sv and neoloop count data
sv_neoloop = pd.read_csv('outputs/sv_neoloop_count.tsv', sep='\t', index_col=0)
sv_neoloop['WGS_ID'] = pd.to_numeric(sv_neoloop['WGS_ID'], errors='coerce').astype('Int64')
sv_neoloop['WGS_ID'] = pd.to_numeric(sv_neoloop['WGS_ID'], errors='coerce').astype('str')
sv_neoloop['HiC_ID'] = pd.to_numeric(sv_neoloop['HiC_ID'], errors='coerce').astype('Int64')
sv_neoloop['HiC_ID'] = pd.to_numeric(sv_neoloop['HiC_ID'], errors='coerce').astype('str')

# merge two datasets
merged = pd.merge(metadata, sv_neoloop, on=['WGS_ID', 'HiC_ID'], how='left')

# prepare lists to make final dictionary for plotting
pr_sv = merged[merged['Tissue Type'] == 'Primary']['SV_count'].tolist()
pr_neo = merged[merged['Tissue Type'] == 'Primary']['neoloop_count'].tolist()
rc_sv = merged[merged['Tissue Type'] == 'Local Regional Recurrence']['SV_count'].tolist()
rc_neo = merged[merged['Tissue Type'] == 'Local Regional Recurrence']['neoloop_count'].tolist()

# create dictionary format for plotting
mydict = {}
mydict['Primary SV'] = pr_sv
mydict['Primary Neoloop'] = pr_neo
mydict['Recurrent SV'] = rc_sv
mydict['Recurrent Neoloop'] = rc_neo

########################################################################################
fontsize = 16

fig, ax = plt.subplots(figsize=(13, 10))

box_colors = ['#DD1717', '#0F4392', '#DD1717', '#0F4392']  # specify colors for boxplots
scatter_colors = ['#DD1717', '#0F4392', '#DD1717', '#0F4392']  # specify colors for scatter plots

data = list(mydict.values())
positions = np.arange(1, len(mydict) + 1) + 0.1  # adjust the offset value here

boxplot = ax.boxplot(data, positions=positions, widths=0.4, patch_artist=True)

# customize box properties
for i, box in enumerate(boxplot['boxes']):
    box.set(facecolor=box_colors[i], alpha=0.5, linewidth=2)

# customize whisker properties
for whisker in boxplot['whiskers']:
    whisker.set(color='gray', linewidth=2, linestyle='--')

# customize median properties
for median in boxplot['medians']:
    median.set(color='black', linewidth=2)

ax.set_xlim(0.5, 4.5)
    
# adjust ytick label size and format
ax.set_yticklabels([int(label) for label in ax.get_yticks()], fontsize=fontsize-2)

# adjust xtick labels and positions
ax.set_xticks([1.5, 3.5])
ax.set_xticklabels(['Primary', 'Recurrent'], fontsize=fontsize)

ax.set_title('Structural Variants and Neoloops in Primary/Recurrent GBM Samples', fontsize=fontsize+2)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# add scatter plots
scatter_offset = 0.25
for i, (key, values) in enumerate(mydict.items()):
    x = np.random.normal(i + 1 - scatter_offset, 0.04, size=len(values))
    ax.scatter(x, values, color=scatter_colors[i], alpha=0.5)

sv_patch = mpatches.Patch(facecolor=box_colors[0], edgecolor='black', alpha=0.5, label='SV Burden')
neoloop_patch = mpatches.Patch(facecolor=box_colors[1], edgecolor='black', alpha=0.5, label='Neoloop Count')
ax.legend(handles=[sv_patch, neoloop_patch], loc='upper right', fontsize=fontsize-2)
    
fig.set_facecolor('white')
plt.show()
