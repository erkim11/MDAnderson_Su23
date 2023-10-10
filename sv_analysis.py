import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import os

# import and clean metadata
metadata = pd.read_excel('CATALYST - MS96 Selected.v5-ATAC.Seq.xlsx')
metadata = metadata[['HiC Seq Corrected', 'WGS Seq ID', 'ecDNA']]
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('Int64')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('str')
metadata['WGS Seq ID'] = pd.to_numeric(metadata['WGS Seq ID'], errors='coerce').astype('Int64')
metadata['WGS Seq ID'] = pd.to_numeric(metadata['WGS Seq ID'], errors='coerce').astype('str')

df = pd.DataFrame()
df['HiC_ID'] = metadata['HiC Seq Corrected']
df['WGS_ID'] = metadata['WGS Seq ID']

########################################################################################
sv_array = np.zeros((8,df.shape[0]), dtype=int)
svs = ['translocation', 'inversion', 'duplication', 'deletion']

# recode sv types to familiar terms
replacing_rules = {
    'EVENTTYPE=BND': 'translocation',
    'EVENTTYPE=DEL': 'deletion',
    'EVENTTYPE=INV': 'inversion',
    'EVENTTYPE=DUP': 'duplication',
    'EVENTTYPE=INS': 'insertion'
}

# parse HiC sv files
for i in range(0, df.shape[0]): # loop through samples
    if df.iloc[i,0] != '<NA>':
        file_path = 'data/'+df.iloc[i,0]+'-HC01.SV.txt'
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0: # make sure file exists and SV exists
            mylist = pd.read_csv(file_path,sep='\t',header=None)[5].tolist()
            for j in range(0, 4): # loop through svs
                sv_count = mylist.count(svs[j])
                sv_array[j, i] = sv_count

# parse WGS sv files
for i in range(0, df.shape[0]): # loop through samples
    if df.iloc[i,1] != '<NA>':
        file_path = 'data_wgs/'+df.iloc[i,1]+'-WG01.sv.vcf'
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0: # make sure file exists and SV exists
            file = pd.read_csv(file_path, sep='\t', comment="#", header=None)
            mylist = []
            for row in range(0, file.shape[0]):
                eventtype = [item for item in file.iloc[row, 7].split(';') if 'EVENTTYPE=' in item]
                mylist.append(replacing_rules.get(eventtype[0], eventtype[0])) # recode using replacing_rules
            for j in range(4, 8): # loop through svs
                sv_count = mylist.count(svs[j-4])
                sv_array[j, i] = sv_count

########################################################################################
ins = np.zeros((1,df.shape[0]), dtype=int)

for i in range(0, df.shape[0]): # loop through samples
    if df.iloc[i,1] != '<NA>':
        file_path = 'data_wgs/'+df.iloc[i,1]+'-WG01.sv.vcf'
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            file = pd.read_csv(file_path, sep='\t', comment="#", header=None)
            mylist = []
            for row in range(0, file.shape[0]):
                eventtype = [item for item in file.iloc[row, 7].split(';') if 'EVENTTYPE=' in item]
                mylist.append(replacing_rules.get(eventtype[0], eventtype[0]))
                sv_count = mylist.count('insertion')
                ins[0, i] = sv_count

########################################################################################
df['HiC_translocation'] = sv_array[0,:]
df['HiC_inversion'] = sv_array[1,:]
df['HiC_duplication'] = sv_array[2,:]
df['HiC_deletion'] = sv_array[3,:]

df['WGS_translocation'] = sv_array[4,:] - ins[0,:] - sv_array[5,:] - sv_array[6,:] - sv_array[7,:]
df['WGS_inversion'] = sv_array[5,:]
df['WGS_duplication'] = sv_array[6,:]
df['WGS_deletion'] = sv_array[7,:]

df['HiC_count'] = df.iloc[:, 2:6].sum(axis=1)
df['WGS_count'] = df.iloc[:, 6:10].sum(axis=1)

df_sorted = df.sort_values(by='WGS_count', ascending=True)
df_sorted = df_sorted[~((df_sorted['HiC_ID'] == '<NA>') | (df['WGS_ID'] == '<NA>'))]
df_sorted = df_sorted.drop(df_sorted.index[0])
df_sorted['ecDNA'] = metadata['ecDNA']

# create dictionary to make plots
hic_dict = {}
wgs_dict = {}

hic_dict['translocation'] = df_sorted['HiC_translocation'].tolist()
hic_dict['inversion'] = df_sorted['HiC_inversion'].tolist()
hic_dict['duplication'] = df_sorted['HiC_duplication'].tolist()
hic_dict['deletion'] = df_sorted['HiC_deletion'].tolist()

wgs_dict['translocation'] = df_sorted['WGS_translocation'].tolist()
wgs_dict['inversion'] = df_sorted['WGS_inversion'].tolist()
wgs_dict['duplication'] = df_sorted['WGS_duplication'].tolist()
wgs_dict['deletion'] = df_sorted['WGS_deletion'].tolist()

########################################################################################
# plot bar plot with sv info from wgs and hic
width = 0.5
fontsize = 32
fig, ax = plt.subplots(nrows=2, figsize=(28, 24), gridspec_kw={'height_ratios': [3, 1]})

colors = ['#5B9EA6', '#A9D4D9', '#D9A384', '#BF6B63']

# Plot for the first subplot
bottom = np.zeros(df_sorted.shape[0])
for i, (num, sv_counts) in enumerate(wgs_dict.items()):
    p = ax[0].bar(df_sorted['WGS_ID'].tolist(), sv_counts, width, label=num, bottom=bottom, color=colors[i])
    bottom += sv_counts

ax[0].set_title("Number of Structural Variants in GBM Samples", fontsize=fontsize+4)

handles, labels = ax[0].get_legend_handles_labels()
labels[0] = 'translocation / breakend'
ax[0].legend(handles, labels, loc="upper left", fontsize=fontsize)
ax[0].set_xticks([])
ax[0].tick_params(axis='y', labelsize=fontsize)

ax[0].text(-1.8, 3350, "Whole Genome Sequencing", fontsize=fontsize)
# ax[0].text(-1.8, 2800, "Whole Genome Sequencing", fontsize=fontsize)

ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)

# samples with ecDNA
for i in range(0, df_sorted.shape[0]):
    if str(df_sorted.iloc[i, 12]) != 'nan':
        ax[0].plot(i-0.05, df_sorted.iloc[i, 11]+100, marker='x', markersize=fontsize/2, color='red')

# plot for the second subplot with reversed bars aligned to the top
bottom = np.zeros(df_sorted.shape[0])
    
for i, (num, sv_counts) in enumerate(hic_dict.items()):
    bottom -= sv_counts
    p = ax[1].bar(df_sorted['HiC_ID'].tolist(), sv_counts, width, label=num, bottom=bottom, color=colors[i])

# ax[1].set_title("Number of Structural Variants in GBM Samples", fontsize=fontsize+2)
ax[1].set_xticks([])
ax[1].tick_params(axis='y', labelsize=fontsize)

ax[1].text(-1.8, -115, "HiC Sequencing", fontsize=fontsize)
# ax[1].text(-1.8, -110, "HiC Sequencing", fontsize=fontsize)

yticks = range(0, -161, -20)  # Set the desired tick positions
yticklabels = [str(abs(ytick)) for ytick in yticks]  # Convert the tick positions to positive values as labels

ax[1].set_yticks(yticks)  # Set the tick positions
ax[1].set_yticklabels(yticklabels)  # Set the tick labels

ax[1].spines['bottom'].set_visible(False)
ax[1].spines['right'].set_visible(False)

# mark samples with ecDNA
for i in range(0, df_sorted.shape[0]):
    if str(df_sorted.iloc[i, 12]) != 'nan':
        ax[1].plot(i-0.05, df_sorted.iloc[i, 10]*-1-7.5, marker='x', markersize=fontsize/2, color='red')

fig.set_facecolor('white')

plt.tight_layout()
plt.show()
