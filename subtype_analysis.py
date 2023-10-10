import numpy as np
import pandas as pd
from biomart import BiomartServer
from scipy import stats
import matplotlib.pyplot as plt
import os

# function for converting gene symbol to ensembl ID
def convert_gene_symbol_to_ensembl(gene_symbols):
    # connect to the Biomart server
    server = BiomartServer("http://www.ensembl.org/biomart")

    # select the Ensembl genes dataset
    ensembl_dataset = server.datasets['hsapiens_gene_ensembl']

    # define the attributes to retrieve
    attributes = ['ensembl_gene_id', 'external_gene_name']

    # perform the query
    response = ensembl_dataset.search({'attributes': attributes, 'filters': {'external_gene_name': gene_symbols}})

    # extract the results
    results = response.content.decode().strip().split('\n')

    # parse the results into a dictionary of gene symbols to ensembl IDs
    symbol_to_ensembl = {}
    for line in results:
        ensembl_id, symbol = line.split('\t')
        symbol_to_ensembl[symbol] = ensembl_id

    return symbol_to_ensembl

# grab signature gene list from centroids data of Verhaak et al.
sig_genes = pd.read_excel('ClaNC840_centroids.xls', header=None)
sig_genes.columns = sig_genes.iloc[2]
sig_genes = sig_genes[3:].reset_index(drop=True)
sig_genes_list = sig_genes.iloc[:,0].tolist()
sig_ensembl = convert_gene_symbol_to_ensembl(sig_genes_list)

# add ensembl id to the dataframe using the new sig_ensembl list
sig_genes['id'] = sig_genes.iloc[:,0].map(sig_ensembl)

########################################################################################
# import and clean metadata
metadata = pd.read_excel('CATALYST - MS96 Selected.v5-ATAC.Seq.xlsx')
metadata = metadata[['HiC Seq Corrected', 'RNA Seq ID']]
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('Int64')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('str')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('Int64')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('str')
metadata = metadata[~metadata['RNA Seq ID'].str.contains('<NA>')]

# list to store subtypes
subtypes = []

for row in range(0, metadata.shape[0]):
    sample_exp_df = pd.DataFrame(list(sig_genes['id']))
    
    rna_sample = metadata.iloc[row, 1]
    new_sample = pd.read_csv('data_rna/'+rna_sample+'-RN01.ReadsPerGene.out.tab', sep='\t', header=None)
    new_sample = new_sample.iloc[4:,:2]
    new_sample.iloc[:, 0] = new_sample.iloc[:, 0].str.split('.').str[0]
    
    merged_df = pd.merge(sample_exp_df, new_sample, on=0, how='left')
    clean_final = merged_df[1]
    clean_final = np.array(clean_final)

    centroid_labels = sig_genes.iloc[:, 2:6]
    distances = np.nansum((clean_final[:, np.newaxis] - centroid_labels.values)**2, axis=0)
    closest_centroids = np.argmin(distances)
    assigned_clusters = ['Proneural', 'Neural', 'Classical', 'Mesenchymal'][closest_centroids]
    
    subtypes.append(assigned_clusters)

########################################################################################
# load save subtype dataset
subtypes_data = pd.read_csv('outputs/subtypes_most_accurate_final.txt', sep='\t', index_col=0)
subtypes_data = subtypes_data.rename(columns={'sample': 'RNA Seq ID'})
subtypes_data = subtypes_data.astype(str)

metadata = pd.merge(metadata, subtypes_data, on='RNA Seq ID')

# add neoloop count info
neoloops = np.zeros(metadata.shape[0], dtype='int64')

for row in range(0, metadata.shape[0]):
    sample = metadata.iloc[row, 0]
    if sample != '<NA>':
        file_path = 'data/'+str(sample)+'-HC01.neo-loops.txt'
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            file = pd.read_csv(file_path, sep='\t', header=None)
            num_neoloop = file.shape[0]
            neoloops[row] = num_neoloop

metadata['Neoloops'] = neoloops

########################################################################################
# count num of neoloops by subtype
metadata_pn = metadata[metadata['subtypes'] == 'Proneural']
pn_loops = metadata_pn['Neoloops'].tolist()
print('pn length: '+str(len(pn_loops)))

metadata_cl = metadata[metadata['subtypes'] == 'Classical']
cl_loops = metadata_cl['Neoloops'].tolist()
print('cl length: '+str(len(cl_loops)))

metadata_mc = metadata[metadata['subtypes'] == 'Mesenchymal']
mc_loops = metadata_mc['Neoloops'].tolist()
print('mc length: '+str(len(mc_loops)))

metadata_nr = metadata[metadata['subtypes'] == 'Neural']
nr_loops = metadata_nr['Neoloops'].tolist()
print('nr length: '+str(len(nr_loops)))

# perform ANOVA test
f_value, p_value = stats.f_oneway(mc_loops, pn_loops, cl_loops)

# create dict for visualization
mydict = {}
mydict['Proneural'] = pn_loops
mydict['Mesenchymal'] = mc_loops
mydict['Classical'] = cl_loops

########################################################################################
# grouped boxplot
fontsize = 16

fig, ax = plt.subplots(figsize=(13, 10))

box_colors = ['#353535', '#D9D9D9', '#3C6E71']
scatter_colors = ['#353535', '#D9D9D9', '#3C6E71']

data = list(mydict.values())
positions = np.arange(1, len(mydict) + 1) + 0.1  # Adjust the offset value here

boxplot = ax.boxplot(data, positions=positions, widths=0.32, patch_artist=True)

# customize box properties
for i, box in enumerate(boxplot['boxes']):
    box.set(facecolor=box_colors[i], alpha=0.7, linewidth=2)

# customize whisker properties
for whisker in boxplot['whiskers']:
    whisker.set(color='gray', linewidth=2, linestyle='--')

# customize median properties
for median in boxplot['medians']:
    median.set(color='black', linewidth=2)

ax.set_xlim(0.5, 3.5)
ax.set_yticklabels([int(label) for label in ax.get_yticks()], fontsize=fontsize)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['Proneural', 'Mesenchymal', 'Classical'],fontsize=fontsize)
ax.set_title('Number of Neoloops per GBM Subtypes', fontsize=fontsize+2)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.text(0.7, 4000, "p-value: "+str(round(p_value, 4)), fontsize=fontsize) 

# add scatter plots
scatter_offset = 0.2
for i, (key, values) in enumerate(mydict.items()):
    x = np.random.normal(i + 1 - scatter_offset, 0.04, size=len(values))
    ax.scatter(x, values, color=scatter_colors[i], alpha=0.5)

fig.set_facecolor('white')
plt.show()

########################################################################################
# prepare data for pie chart
loop_count = []
loop_count.append(np.sum(pn_loops))
loop_count.append(np.sum(mc_loops))
loop_count.append(np.sum(cl_loops))

labels = ['Proneural', 'Mesenchymal', 'Classical']
colors = ['#353535', '#D9D9D9', '#3C6E71']

# plot pie chart
fig, ax = plt.subplots(figsize=(10, 10))

patches, texts, autotexts = ax.pie(loop_count, labels=labels, autopct='%1.1f%%',
                                   colors = colors, textprops={'fontsize': fontsize})

for patch in patches:
    patch.set_facecolor(color=(patch.get_facecolor()[:-1] + (0.8,)))

ax.set_title('Distribution of Neoloops Across GBM Subtypes', fontsize=fontsize+2)
fig.set_facecolor('white')

plt.show()
