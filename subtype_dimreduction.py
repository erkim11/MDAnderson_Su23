import numpy as np
import pandas as pd
import umap
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from biomart import BiomartServer
import matplotlib.pyplot as plt
import os

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

    # parse the results into a list of Ensembl IDsd
    symbol_to_ensembl = []
    for line in results:
        ensembl_id, symbol = line.split('\t')
        symbol_to_ensembl.append(ensembl_id)

    return symbol_to_ensembl

sig_genes = pd.read_excel('ClaNC840_centroids.xls')
sig_genes_list = sig_genes.iloc[2:,0].tolist()

ensembl_list = convert_gene_symbol_to_ensembl(sig_genes_list)

metadata = pd.read_excel('CATALYST - MS96 Selected.v5-ATAC.Seq.xlsx')
metadata = metadata[['HiC Seq Corrected', 'RNA Seq ID']]

metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('Int64')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('str')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('Int64')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('str')

########################################################################################
# import rna data
rnadata = pd.read_csv('data_rna/235921-RN01.ReadsPerGene.out.tab', sep='\t', header=None)
rnadata = rnadata.iloc[4:]
rnadata.iloc[:, 0] = rnadata.iloc[:, 0].str.split('.').str[0]
filtered_rnadata = rnadata[rnadata.iloc[:, 0].isin(ensembl_list)]
filtered_rnadata.iloc[:, 1]

mydata = np.zeros((metadata.shape[0], 737), dtype=int)

for row in range(0, metadata.shape[0]):
    rna_sample = metadata.iloc[row, 1]
    rna_data = pd.read_csv('data_rna/'+rna_sample+'-RN01.ReadsPerGene.out.tab', sep='\t', header=None)
    
    rna_data = rna_data.iloc[4:] # remove first four rows
    
    rna_data.iloc[:, 0] = rna_data.iloc[:, 0].str.split('.').str[0]
    filtered_rna_data = rna_data[rna_data.iloc[:, 0].isin(ensembl_list)]
    final_rna_data = filtered_rna_data.iloc[:, 1]
    
    rna_list = final_rna_data.tolist()
    for gene in range(0, 737):
        mydata[row][gene] = rna_list[gene]

column_sums = mydf.sum(axis=0)
sum_row = pd.DataFrame(column_sums).T
sum_row.index = ['sum']
mydf_with_sum = pd.concat([mydf, sum_row])

sorted_sum_row = mydf_with_sum.loc['sum'].sort_values(ascending=False)
top_500_values = sorted_sum_row.head(600)
top_500_columns = top_500_values.index
filtered_df = mydf_with_sum[top_500_columns]

final_df = filtered_df.drop(filtered_df.index[-1])

########################################################################################
# perform t-SNE analysis
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(final_df)

# plot the t-SNE results
fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(tsne_result[:, 0], tsne_result[:, 1])
ax.set_xlabel('t-SNE Dimension 1')
ax.set_ylabel('t-SNE Dimension 2')
ax.set_title('t-SNE Analysis')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.set_facecolor('white')

plt.show()

########################################################################################
# create a PCA object and specify the number of components
pca = PCA(n_components=2)

# fit the PCA model to the data
pca.fit(final_df)

# transform the data to the principal components
X_pca = pca.transform(final_df)

# access the explained variance ratio
explained_variance_ratio = pca.explained_variance_ratio_

# plot the transformed data
fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(X_pca[:, 0], X_pca[:, 1])
ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.set_title('PCA Analysis')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.set_facecolor('white')
plt.show()

########################################################################################
# perform UMAP analysis
umap_model = umap.UMAP(n_components=2, random_state=42)

umap_result = umap_model.fit_transform(final_df)

fig, ax = plt.subplots(figsize=(10, 8))

ax.scatter(umap_result[:, 0], umap_result[:, 1])
ax.set_xlabel('UMAP Dimension 1')
ax.set_ylabel('UMAP Dimension 2')
ax.set_title('UMAP Analysis')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.set_facecolor('white')
plt.show()
