import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.preprocessing import RobustScaler
import os

# load 840 signature genes data
sig_genes = pd.read_csv('outputs/sig_genes_Verhaak.txt', sep='\t', index_col=0)

sig_genes['Subtype'] = sig_genes[['Proneural', 'Classical', 'Mesenchymal']].idxmax(axis=1)
sig_genes = sig_genes.sort_values(by='Subtype', ascending=False)

# load subtype data
subtypes = pd.read_csv('outputs/subtypes_most_accurate_final.txt', sep='\t', index_col=0)
subtypes = subtypes.sort_values(by='subtypes', ascending=False)

########################################################################################
# create dataframe of gbm samples and gene expression levels
mydata = np.zeros((subtypes.shape[0], sig_genes.shape[0]))

for i in range(0, subtypes.shape[0]): # i is sample index
    sample = subtypes.iloc[i, 1]
    rna_data = pd.read_csv('data_rna/'+str(sample)+'-RN01.ReadsPerGene.out.tab', sep='\t', header=None)
    rna_data.iloc[:, 0] = rna_data.iloc[:, 0].str.split('.').str[0]
    
    for j in range(0, sig_genes.shape[0]): # j is gene index
        gene_id = sig_genes.iloc[j, 6]
        if len(rna_data[rna_data.iloc[:, 0] == gene_id]) != 0:
            mydata[i, j] = rna_data[rna_data.iloc[:, 0] == gene_id][1].values[0]

row_indices = subtypes.iloc[:, 1].tolist()
column_indices = sig_genes.iloc[:, 0].tolist()

mydf = pd.DataFrame(mydata, columns=column_indices, index=row_indices)

########################################################################################
# sort by gene
pr_column_sums = mydf.iloc[0:35, 0:450].sum(axis=0, skipna=True)
pr_sum_row = pd.DataFrame(pr_column_sums).T
pr_sum_row.index = ['pr_sum']

ms_column_sums = mydf.iloc[35:59, 450:614].sum(axis=0, skipna=True)
ms_sum_row = pd.DataFrame(ms_column_sums).T
ms_sum_row.index = ['ms_sum']

cl_column_sums = mydf.iloc[59:79, 614:840].sum(axis=0, skipna=True)
cl_sum_row = pd.DataFrame(cl_column_sums).T
cl_sum_row.index = ['cl_sum']

mydf = pd.concat([mydf, pr_sum_row])
mydf = pd.concat([mydf, ms_sum_row])
mydf = pd.concat([mydf, cl_sum_row])

sorted_df1 = mydf.iloc[:, :450].T.sort_values(by='pr_sum', ascending=False).T
sorted_df2 = mydf.iloc[:, 450:614].T.sort_values(by='ms_sum', ascending=False).T
sorted_df4 = mydf.iloc[:, 614:840].T.sort_values(by='cl_sum', ascending=False).T

concatenated_df = pd.concat([sorted_df1, sorted_df2, sorted_df4], axis=1)

########################################################################################
# sort by sample
pr_row_sums = mydf.iloc[0:35, 0:450].sum(axis=1, skipna=True)
pr_sum_col = pd.DataFrame(pr_row_sums)

ms_row_sums = mydf.iloc[35:59, 450:614].sum(axis=1, skipna=True)
ms_sum_col = pd.DataFrame(ms_row_sums)

cl_row_sums = mydf.iloc[59:79, 614:840].sum(axis=1, skipna=True)
cl_sum_col = pd.DataFrame(cl_row_sums)

mydf['pr_sum'] = pr_sum_col
mydf['ms_sum'] = ms_sum_col
mydf['cl_sum'] = cl_sum_col

sorted_df1 = mydf.iloc[0:35, :].sort_values(by='pr_sum', ascending=False)
sorted_df2 = mydf.iloc[35:59, :].sort_values(by='ms_sum', ascending=False)
sorted_df4 = mydf.iloc[59:79, :].sort_values(by='cl_sum', ascending=False)

concatenated_df = pd.concat([sorted_df1, sorted_df2, sorted_df4])

# cut out the sum columns and rows
concatenated_df = concatenated_df.iloc[:79,:840]

# normalize data
final_array = concatenated_df.iloc[0:79,:].values
normalized_array = preprocessing.scale(final_array)

########################################################################################
# define heatmap function
def heatmap(data, ax=None, **kwargs):
    '''
    Create a heatmap from a numpy array and two lists of labels.
    '''

    if ax is None:
        ax = plt.gca()

    # rotate the data 90 degrees
    data = np.transpose(data)

    im = ax.imshow(data, **kwargs)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

    ax.set_aspect(0.05)
    
    ax.axvline(34.5, color='black', linewidth=1)
    ax.axvline(58.5, color='black', linewidth=1)
    
    ax.set_title('Gene Expression Levels of GBM Samples', fontsize=20, pad=45)
    ax.text(14.5, -15, 'Proneural', fontsize=18)
    ax.text(42.5, -15, 'Mesenchymal', fontsize=18)
    ax.text(66, -15, 'Classical', fontsize=18)
    
    return im

# plot heatmap
fig, ax = plt.subplots(figsize=(20, 20))
fig.set_facecolor('white')

im = heatmap(data=normalized_array, cmap='RdBu_r', vmin=np.min(normalized_array)+0.2, vmax=np.max(normalized_array)-5.8)
