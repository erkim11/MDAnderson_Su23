import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import pyth.plugins.rtf15.reader as reader
import os

# import metadata
metadata = pd.read_excel('CATALYST - MS96 Selected.v5-ATAC.Seq.xlsx')
metadata = metadata[['HiC Seq Corrected', 'WGS Seq ID', 'ecDNA']]

metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('Int64')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('str')
metadata['WGS Seq ID'] = pd.to_numeric(metadata['WGS Seq ID'], errors='coerce').astype('Int64')
metadata['WGS Seq ID'] = pd.to_numeric(metadata['WGS Seq ID'], errors='coerce').astype('str')

# import IDH mutant status data
idh = pd.read_csv('idh.txt', sep='\t', header=None)
idh = idh.rename(columns={0: 'WGS_ID'})
idh['WGS_ID'] = pd.to_numeric(idh['WGS_ID'], errors='coerce').astype('str')

df = pd.DataFrame()
df['HiC_ID'] = metadata['HiC Seq Corrected']
df['WGS_ID'] = metadata['WGS Seq ID']

########################################################################################
# merge dataframes
mydf = pd.merge(df, idh, how='left', on='WGS_ID')

# add column for number of neoloops in each sample
neoloops = np.zeros(mydf.shape[0], dtype='int64')

for row in range(0, mydf.shape[0]):
    sample = mydf.iloc[row, 0]
    if sample != '<NA>':
        file_path = 'data/'+sample+'-HC01.neo-loops.txt'
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            file = pd.read_csv(file_path, sep='\t', header=None)
            num_neoloop = file.shape[0]
            neoloops[row] = num_neoloop

mydf['neoloops'] = neoloops # add neoloops list as column in mydf

mydf = mydf.dropna()
mydf = mydf[~mydf['HiC_ID'].str.contains('<NA>')]

# create mut_loops list that contain the number of neoloops in IDH mutant samples
mydf_mut = mydf[mydf[1] == 'mut']
mut_loops = mydf_mut['neoloops'].tolist()

# create wt_loops list that contain the number of neoloops in IDH wild-type samples
mydf_wt = mydf[mydf[1] == 'wt']
wt_loops = mydf_wt['neoloops'].tolist()

# make it into a dictionary format for plotting
mydict = {}
mydict['mut'] = mut_loops
mydict['wt'] = wt_loops

########################################################################################
fontsize = 16

# plot boxplot with scatter
fig, ax = plt.subplots(figsize=(15, 12))

box_colors = ['#DD1717', '#0F4392']  # specify colors for boxplots
scatter_colors = ['#DD1717', '#0F4392']  # specify colors for scatter plots

data = list(mydict.values())
positions = np.arange(1, len(mydict) + 1) + 0.1  # adjust the offset value here

boxplot = ax.boxplot(data, positions=positions, widths=0.32, patch_artist=True)

# customize box properties
for i, box in enumerate(boxplot['boxes']):
    box.set(facecolor=box_colors[i], alpha=0.5, linewidth=2)

# customize whisker properties
for whisker in boxplot['whiskers']:
    whisker.set(color='gray', linewidth=2, linestyle='--')

# customize median properties
for median in boxplot['medians']:
    median.set(color='black', linewidth=2)

ax.set_xlim(0.5, 2.5)
    
# adjust ytick label size and format
ax.set_yticklabels([int(label) for label in ax.get_yticks()], fontsize=fontsize)

# adjust xtick labels and positions
ax.set_xticks([1, 2])
ax.set_xticklabels(['IDH Mutant', 'IDH Wild-Type'], fontsize=fontsize)

ax.set_title('Number of Neoloops in GBM Samples', fontsize=fontsize+2)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.text(0.62, 4000, "p-value: "+str(round(p_value, 4)), fontsize=fontsize) 

# scatter plots
scatter_offset = 0.2 # offset parameter
for i, (key, values) in enumerate(mydict.items()):
    x = np.random.normal(i + 1 - scatter_offset, 0.04, size=len(values))
    ax.scatter(x, values, color=scatter_colors[i], alpha=0.5)

fig.set_facecolor('white')
plt.show()
