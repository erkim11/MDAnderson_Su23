import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# list of sample hic id
samples = ['236433', '236434', '236435', '236436', '236437', '236438', '236439',
           '236440', '236441', '236442', '236444', '236452', '236453',
           '236454', '236455', '236456', '236457', '236458', '236459', '236460',
           '236461', '236462', '236463',
           '231038', '231039', '231040', '231041', '231042', '231043', '231044',
           '231045', '231046', '231047', '231048', '231049', '231050', '231051',
           '231401', '231402', '231403', '231404', 
           '231407', '231408', '231409', '231410', '235192', '235193',
           '235195', '235197', '235198', '235199', '235200',
           '235202', '235203', '235204', '235205', '235206', '235207', '235208',
           '235209', '235210', '235211', '235212',
           '236445', '236447', '236448', '236449', '236450', '236451', '236464',
           '236465', '236466', '236467', '236496', '236497', '236498', '236499',
           '236500', '236501']

# final result dictionary
overlapping_dict = {} # key is sample id and value is a list of overlapping genes
flanking_dict = {} # key is sample id and value is a list of flanking genes

flanking_distance = 100000 # parameter open to adjust

# load list of oncogenes from cancer gene census
oncogenes_file = pd.read_csv('data/CancerGeneCensus.tsv',sep='\t',header=None)
oncogenes_df = pd.DataFrame(oncogenes_file.iloc[1:, [0]])
oncogenes_df = oncogenes_df.rename(columns={0: 'gene'})

oncogenes_df['chr'] = oncogenes_file.iloc[:, 3].str.split(':', expand=True)[0]
oncogenes_df['start'] = oncogenes_file.iloc[:, 3].str.split(':', expand=True)[1].str.split('-', expand=True)[0]
oncogenes_df['end'] = oncogenes_file.iloc[:, 3].str.split(':', expand=True)[1].str.split('-', expand=True)[1]

for sample in samples:
    neoloops_file = pd.read_csv('data/'+sample+'-HC01.neo-loops.txt',sep='\t',header=None) # load neoloop file
    
    # tidy neoloop data
    neoloops_df = pd.DataFrame(neoloops_file.iloc[0:, [0]])
    neoloops_df = neoloops_df.rename(columns={0: 'chr_1'})
    neoloops_df['chr_1'] = neoloops_df['chr_1'].str.replace('chr', '')
    neoloops_df['start_1'] = neoloops_file.iloc[:, 1]
    neoloops_df['end_1'] = neoloops_file.iloc[:, 2]
    neoloops_df['chr_2'] = neoloops_file.iloc[:, 3].str.replace('chr', '')
    neoloops_df['start_2'] = neoloops_file.iloc[:, 4]
    neoloops_df['end_2'] = neoloops_file.iloc[:, 5]

# make two lists to store name of genes where neoloops are overlapping and/or flanking
    overlapping_genes = []
    flanking_genes = []
    
    for i in range(0, oncogenes_df.shape[0]):
        for j in range(0, neoloops_df.shape[0]):
            neoloop_chr = neoloops_df.iat[j, 0]
            neoloop_start = int(neoloops_df.iat[j, 1])
            neoloop_end = int(neoloops_df.iat[j, 5])

            gene_name = oncogenes_df.iat[i, 0]
            gene_chr = oncogenes_df.iat[i, 1]
            gene_start = oncogenes_df.iat[i, 2]
            gene_end = oncogenes_df.iat[i, 3]

            if gene_start != '' and gene_end != '':
                gene_start = int(gene_start)
                gene_end = int(gene_end)

                if neoloop_chr == gene_chr:
                    if (gene_start - neoloop_end)*(gene_end - neoloop_start) < 0:
                        overlapping_genes.append(gene_name)
                    elif gene_end < neoloop_start and gene_end > (neoloop_start - flanking_distance):
                        flanking_genes.append(gene_name)
                    elif gene_start < (neoloop_end + flanking_distance) and gene_start > neoloop_end:
                        flanking_genes.append(gene_name)
                        
    overlapping_dict[sample] = overlapping_genes
    flanking_dict[sample] = flanking_genes
    
    print(sample+'-HC01 analyze successful')
 
max_length1 = max(len(arr) for arr in overlapping_dict.values())
max_length2 = max(len(arr) for arr in flanking_dict.values())

# pad the shorter arrays with None to make them equal in length
for key, arr in overlapping_dict.items():
    if len(arr) < max_length1:
        overlapping_dict[key] += [None] * (max_length1 - len(arr))

for key, arr in flanking_dict.items():
    if len(arr) < max_length2:
        flanking_dict[key] += [None] * (max_length2 - len(arr))

# save the DataFrame as a text file
df1 = pd.DataFrame.from_dict(overlapping_dict)
df2 = pd.DataFrame.from_dict(flanking_dict)

df1.to_csv('overlapping.txt', sep='\t', index=False)
df2.to_csv('flanking.txt', sep='\t', index=False)
