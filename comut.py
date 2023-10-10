import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from comut import comut
from comut import fileparsers

########################################################################################
# version with the indicators for same patients

data = pd.read_csv('outputs/comut_data2.csv', sep=',')
# tidy the data for comut plot
data = pd.melt(data, id_vars=['sample', 'group', 'svs'], value_vars=['Diagnosis', 'Subtype', 'IDH', 'ecDNA', 'Tissue', 'Neoloops'],
                  var_name='category', value_name='value')
data['value'].fillna('NA', inplace=True)

# adjust font size
custom_rcParams = {'font.size': 18}
rcParams.update(custom_rcParams)

# color mapping dictionary for categorical data
mapping =  # Diagnosis
          {'Oligodendroglioma': '#D96738',
          'Glioblastoma': '#EBCD51',
          'Astrocytoma': '#48B0C7',
           # Subtype
          'Proneural': '#F0ADA7',
          'Classical': '#6BB2CC',
          'Mesenchymal': '#17517E',
           'NA': 'lightgray',
           # Tissue
          'Primary': '#214785',
          'Recurrent': '#D96767',
           # IDH
          'Mutant': '#E8C9C7',
          'Wild-Type': '#149684',
           # ecDNA
          'Yes': '#473857',
          'No': 'whitesmoke'}

category_order = ['IDH', 'ecDNA', 'Tissue', 'Subtype', 'Diagnosis']
value_range = (0, 700) # range for continuous data

# initialize comut object
toy_comut = comut.CoMut()

# add indicators
indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 5}
toy_comut.add_sample_indicators(data.iloc[:82, :2], name = 'Same patient', plot_kwargs = indicator_kwargs)

# add continuous data
toy_comut.add_continuous_data(data.iloc[410:,:], name = 'Neoloops', 
                              mapping = 'Reds', value_range = value_range)

# add categorical data
toy_comut.add_categorical_data(data.iloc[164:246,:], name = 'IDH', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[246:328, :], name = 'ecDNA', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[328:410,:], name = 'Tissue', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[82:164,:], name = 'Subtype', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[0:82, :], name = 'Diagnosis', mapping = mapping)

# add bar data
toy_comut.add_bar_data(data.iloc[:82,[0, 2]], name = 'SV Burden',
                      mapping = {'svs': '#EE7C58'}, ylabel = 'SV Burden')

# create comut plot
toy_comut.plot_comut(figsize = (20,12), x_padding = 0.04, y_padding = 0.04,
                     hspace = 0.2, heights = {'SV Burden': 4, 'Same patient': 0.5})
toy_comut.figure.set_facecolor('white')
toy_comut.add_unified_legend(bbox_to_anchor = (1, 1.28))
toy_comut.axes['Same patient'].set_xticklabels([])

########################################################################################
# version without the indicators but with color bar instead for same patients

data = pd.read_csv('outputs/comut_data4.csv', sep=',')
# tidy the data
data = pd.melt(data, id_vars=['sample', 'svs'], value_vars=['Diagnosis', 'Subtype', 'IDH', 'ecDNA', 'Tissue', 'Neoloops', 'group'],
                  var_name='category', value_name='value')
data['value'].fillna('NA', inplace=True)
data['category'] = data['category'].replace('group', 'Patient')

# adjust font size
custom_rcParams = {'font.size': 18}
rcParams.update(custom_rcParams)

# color mapping dictionary for categorical data
mapping = {'Oligodendroglioma': '#D96738',
          'Glioblastoma': '#EBCD51',
          'Astrocytoma': '#48B0C7',
           # Subtype
          'Proneural': '#F0ADA7',
          'Classical': '#6BB2CC',
          'Mesenchymal': '#17517E',
           'NA': 'lightgray',
           # Tissue
          'Primary': '#214785',
          'Recurrent': '#D96767',
           # IDH
          'Mutant': '#E8C9C7',
          'Wild-Type': '#149684',
           # ecDNA
          'Yes': '#473857',
          'No': 'whitesmoke'}

category_order = ['IDH', 'ecDNA', 'Tissue', 'Subtype', 'Diagnosis']
value_range = (0, 700)
value_range2 = (0, 100)

# initialize comut object
toy_comut = comut.CoMut()

# add continuous data
toy_comut.add_continuous_data(data.iloc[474:553,:], name = 'Patient', 
                              mapping = 'tab20', value_range = value_range2)
toy_comut.add_continuous_data(data.iloc[395:474,:], name = 'Neoloops', 
                              mapping = 'Reds', value_range = value_range)

# add categorical data
toy_comut.add_categorical_data(data.iloc[158:237,:], name = 'IDH', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[237:316, :], name = 'ecDNA', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[316:395,:], name = 'Tissue', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[79:158,:], name = 'Subtype', mapping = mapping)
toy_comut.add_categorical_data(data.iloc[0:79, :], name = 'Diagnosis', mapping = mapping)

# add bar data
toy_comut.add_bar_data(data.iloc[:79,:2], name = 'SV Burden',
                      mapping = {'svs': '#EE7C58'}, ylabel = 'SV Burden')

toy_comut.plot_comut(figsize = (20,12), x_padding = 0.04, y_padding = 0.04,
                     hspace = 0.2, heights = {'SV Burden': 4, 'Patient': 0.3})
toy_comut.figure.set_facecolor('white')
toy_comut.add_unified_legend(bbox_to_anchor = (1, 0.63))
toy_comut.axes['Patient'].set_xticklabels([])
