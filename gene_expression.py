import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats

genes = ['MLLT6','SBDS','LASP1','ELN','PAX8','CDK12','EGFR',
         'IKZF3','ERBB2','HIP1','IL21R','PML','KRAS','RARA',
         'SMARCE1','CCR7','HOXC13','DDIT3','HOXC11','NCOR1']

overlapping = pd.read_csv('neoloop_gene_overlapping.txt',sep='\t',header=None)
flanking = pd.read_csv('neoloop_gene_flanking.txt',sep='\t',header=None)

########################################################################################
# clean metadata
metadata = pd.read_excel('CATALYST - MS96 Selected.v5-ATAC.Seq.xlsx')
metadata = metadata[['RNA Seq ID', 'HiC Seq Corrected']]
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('Int64')
metadata['RNA Seq ID'] = pd.to_numeric(metadata['RNA Seq ID'], errors='coerce').astype('str')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('Int64')
metadata['HiC Seq Corrected'] = pd.to_numeric(metadata['HiC Seq Corrected'], errors='coerce').astype('str')

########################################################################################
# create gene expression graph for each identified neoloop genes
for gene in genes:
    
    samples_o = [] # samples where select gene is overlapping a neoloop
    samples_f = [] # samples where select gene is flanking a neoloop
    
    for i in range(0, overlapping.shape[1]):
        if gene in overlapping.iloc[:, i].tolist():
            samples_o.append(overlapping.iloc[0, i])

    for i in range(0, flanking.shape[1]):
        if gene in flanking.iloc[:, i].tolist():
            samples_f.append(flanking.iloc[0, i])

# load gene expression read count dataset
    gene_reads = pd.read_csv('data_rna/'+gene+'.txt',sep='\t',header=None)
    gene_reads = gene_reads.drop(2, axis=1)
    gene_reads = gene_reads.drop(3, axis=1)
    gene_reads['RNA Seq ID'] = gene_reads.iloc[:, 0].str.split('-', expand=True)[0]
    gene_reads = gene_reads.drop(0, axis=1)

# add RNA seq id info by left merge
    gene_reads_merged = pd.merge(gene_reads, metadata, on='RNA Seq ID', how='left')
    gene_reads_merged = gene_reads_merged[~gene_reads_merged['HiC Seq Corrected'].str.contains('<NA>')]

# overlapping_status is a list with length same as number of samples
# it has zeros and ones based on whether the sample is within overlapping or not
    overlapping_status = np.zeros(gene_reads_merged.shape[0])
    for row in range(0, gene_reads_merged.shape[0]):
        if gene_reads_merged.iloc[row, 2] in samples_o:
            overlapping_status[row] = True
        else:
            overlapping_status[row] = False

    flanking_status = np.zeros(gene_reads_merged.shape[0])
    for row in range(0, gene_reads_merged.shape[0]):
        if gene_reads_merged.iloc[row, 2] in samples_f:
            flanking_status[row] = True
        else:
            flanking_status[row] = False
    
    gene_reads_merged['Overlapping'] = overlapping_status
    gene_reads_merged['Flanking'] = flanking_status

# make x and y dataset for grouped scatter plot
    x = []
    y = []

    x = x + [0] * len(gene_reads_merged[gene_reads_merged['Overlapping'] == 1][1].tolist())
    y = y + gene_reads_merged[gene_reads_merged['Overlapping'] == 1][1].tolist()

    x = x + [1] * len(gene_reads_merged[gene_reads_merged['Flanking'] == 1][1].tolist())
    y = y + gene_reads_merged[gene_reads_merged['Flanking'] == 1][1].tolist()

    x = x + [2] * len(gene_reads_merged[(gene_reads_merged['Overlapping'] == 0) & (gene_reads_merged['Flanking'] == 0)][1].tolist())
    y = y + gene_reads_merged[(gene_reads_merged['Overlapping'] == 0) & (gene_reads_merged['Flanking'] == 0)][1].tolist()

# prepare sample for t test
    sample1 = gene_reads_merged[gene_reads_merged['Overlapping'] == 1][1].tolist()
    sample2 = gene_reads_merged[gene_reads_merged['Flanking'] == 1][1].tolist()
    sample3 = gene_reads_merged[(gene_reads_merged['Overlapping'] == 0) & (gene_reads_merged['Flanking'] == 0)][1].tolist()

# perform t test
    t_statistic_1, p_value_1 = stats.ttest_ind(sample1, sample3)
    t_statistic_2, p_value_2 = stats.ttest_ind(sample2, sample3)
    
    fontsize = 12

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, color='#173F5F', alpha=0.5)

    ax.set_title(gene+' Expression Based on Neoloop Status', fontsize=fontsize+2)

    # Customize the y-axis tick labels with the original categorical values
    ax.set_xlim(-0.5, 2.5)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Overlapping', 'Flanking', 'Neither'], fontsize=fontsize)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.set_facecolor('white')
    
    ax.text(0.09, -0.14, "p-value: "+str(round(p_value_1, 3))+"                          p-value: "+str(round(p_value_2, 3)), transform=ax.transAxes)
    
    plt.savefig('plots/gene_exp/'+gene+'.png')
