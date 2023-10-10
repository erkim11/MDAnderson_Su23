import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import QuantileTransformer
from matplotlib.colors import ListedColormap
from matplotlib import gridspec
import matplotlib.ticker as mticker
import mpltern
import os

geneexp = pd.read_csv('outputs/subtype_geneexp_corrected_2.txt', sep='\t', index_col=0)
scaler = StandardScaler()

normalized_data = pd.DataFrame(scaler.fit_transform(geneexp), columns=geneexp.columns)

samples_ordered = geneexp.index.tolist()

sig_genes = pd.read_csv('outputs/sig_genes_Verhaak.txt', sep='\t', index_col=0)

########################################################################################
coords = np.zeros((79, 3))
new_subtypes = []

sig_profile = sig_genes.iloc[:, 0:6]

for row in range(0, normalized_data.shape[0]):
    # grab the gene expression profile for each sample
    exp_profile = pd.DataFrame(normalized_data.iloc[row, :]).reset_index()
    exp_profile = exp_profile.rename(columns={'index': 'Gene Symbol'})

    # merge dataframe so that select genes are same
    merged_df = pd.merge(exp_profile, sig_profile, on='Gene Symbol', how='left')
    
    # grab final calculation matrices from merged dataframe
    exp_final = np.array(merged_df.iloc[:, 1])
    centroids = merged_df.iloc[:, [3, 5, 6]]

    # cosine distance
#     distances = distance.cdist(centroids.values.T, exp_final[:, np.newaxis].T, metric='cosine')#     sum_distances = np.sum(distances)
    
    # euclidean distance
    distances = np.nansum((exp_final[:, np.newaxis] - centroids)**2, axis=0)

    coords[row, 0] = distances[0]
    coords[row, 1] = distances[1]
    coords[row, 2] = distances[2]
    
    closest_centroids = np.argmin(distances)
    assigned_clusters = ['Proneural', 'Classical', 'Mesenchymal'][closest_centroids]
    
    new_subtypes.append(assigned_clusters)

new_subtypes = np.array(new_subtypes)

########################################################################################
coords_df = pd.DataFrame(coords, columns=['Proneural', 'Classical', 'Mesenchymal'], index=samples_ordered)

# process coords_df data
coords_df['sum'] = coords_df.sum(axis=1)
coords_df.iloc[:, 0:3] = coords_df.iloc[:, 0:3].div(coords_df['sum'], axis=0)
coords_df = coords_df.iloc[:, :3]
# right now lower value means higher percentage of subtype
# multiply by -1 to make it more intuitive
coords_df = coords_df.multiply(-1)

# scale data
scaler = QuantileTransformer()
normalized_coords = pd.DataFrame(scaler.fit_transform(coords_df), columns=coords_df.columns, index=coords_df.index)

normalized_coords['RNA Seq ID'] = coords_df.index
merged = pd.merge(normalized_coords, subtypes.iloc[:,1:3], on='RNA Seq ID')
subtype_ordered = merged['Subtype']

########################################################################################
mapping = {'Classical': '#3C6E71', 'Mesenchymal': '#D9D9D9', 'Proneural': '#353535'}

# recode the values in subtypes array to colors
new_subtypes_recoded = np.where(new_subtypes == 'Classical', mapping['Classical'],
                np.where(new_subtypes == 'Mesenchymal', mapping['Mesenchymal'],
                np.where(new_subtypes == 'Proneural', mapping['Proneural'], new_subtypes)))

########################################################################################
# prepare ternary coordinates
t0 = normalized_coords['Proneural']
l0 = normalized_coords['Classical']
r0 = normalized_coords['Mesenchymal']

# plot ternary plot
fig = plt.figure(figsize=(8, 8))

ax = fig.add_subplot(projection="ternary")
pc = ax.scatter(t0, l0, r0, color=new_subtypes_recoded)

fontsize=14
ax.grid(linestyle='--')

ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# adjust axis ticks
ax.taxis.set_ticks(ticks, labels=ticks, fontsize=fontsize-2)
ax.laxis.set_ticks(ticks, labels=ticks, fontsize=fontsize-2)
ax.raxis.set_ticks(ticks, labels=ticks, fontsize=fontsize-2)
formatter = mticker.StrMethodFormatter('{x:.0%}')  # format to percentage
ax.taxis.set_major_formatter(formatter)
ax.laxis.set_major_formatter(formatter)
ax.raxis.set_major_formatter(formatter)
ax.tick_params(labelrotation='horizontal')

# adjust axis labels
ax.set_tlabel('Proneural', fontsize=fontsize)
ax.set_llabel('Classical', fontsize=fontsize)
ax.set_rlabel('Mesenchymal', fontsize=fontsize)
ax.taxis.set_label_rotation_mode('horizontal')
ax.laxis.set_label_rotation_mode('horizontal')
ax.raxis.set_label_rotation_mode('horizontal')

# create a legend based on colors
unique_colors = np.unique(new_subtypes_recoded)
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10) for color in unique_colors]
legend_labels = ['Proneural', 'Classical', 'Mesenchymal']  # replace with subtype labels

fig.legend(legend_elements, legend_labels, loc='center right', fontsize=fontsize-2)

plt.title('GBM Subtype Classification', fontsize=fontsize, pad=35)
fig.set_facecolor('white')
plt.show()
