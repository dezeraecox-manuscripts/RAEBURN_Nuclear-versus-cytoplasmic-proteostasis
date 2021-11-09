import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io
import functools

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import OK')

# define location parameters
input_folder = f'results/chaperone_localisation/pixel_collection/'
output_folder = f'results/chaperone_localisation/summary_calculations/'

fret_channel = 3
overlap_threshold = 0.5

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# read in calculated pixel data
file_list = [filename for filename in os.listdir(input_folder) if '.csv' in filename]
pixels = {filename.replace('.csv', ''): pd.read_csv(f'{input_folder}{filename}') for filename in file_list}
pixels.update({key: value.drop([col for col in value.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, value in pixels.items()})

# generate summary df, collect only channel of interest
pixels_compiled = pd.concat(pixels.values())
pixels_compiled = pd.pivot_table(pixels_compiled, index=['x', 'y', 'mask_type', 'cell'], columns=['channel'], values=['intensity']).reset_index()
pixels_compiled.columns = [
    '_'.join(str(val) for val in x) if type(x[1]) == int else x[0]
    for x in pixels_compiled.columns
]


# Add label if aggregate inside unmasked (i.e. same compartment)
aggregate_cells = pixels_compiled[pixels_compiled['mask_type'] == 'aggregate']['cell'].unique()


# generate mean values for each ROI for each timepoint
pixels_mean = pixels_compiled.copy().groupby(['cell', 'mask_type']).median().reset_index()

# assign identifiers
pixels_mean[['treatment', 'chaperone', 'image_number', 'discard2', 'cell_number']] = pixels_mean['cell'].str.split('_', expand=True)
pixels_mean.drop(['discard2'], axis=1, inplace=True)

pixels_mean['aggregate_cell'] = [1 if cell in aggregate_cells else np.nan for cell in pixels_mean['cell']]

# generate nucleus vs cytoplasm ratio
quant_col = 'intensity_3'
info_cols = ['treatment', 'chaperone', 'image_number', 'cell_number', 'aggregate_cell']
cytoplasm = pixels_mean[pixels_mean['mask_type'] == 'cytoplasm'].copy().set_index(info_cols)[quant_col].reset_index().rename(columns={'intensity_3': 'cytoplasm'})
nucleus = pixels_mean[pixels_mean['mask_type'] == 'nucleus'].copy().set_index(info_cols)[quant_col].reset_index().rename(columns={'intensity_3': 'nucleus'})

ratio = functools.reduce(lambda left, right: pd.merge(left, right, on=info_cols, how='outer'), [cytoplasm, nucleus])
ratio['nuc-cyto_ratio'] = ratio['nucleus'] / ratio['cytoplasm']


# save to excel
FileHandling.df_to_excel(output_path=f'{output_folder}summary_calculations.xlsx', sheetnames=['summary', 'nuc-cyto_ratio'], data_frames=[pixels_mean, ratio])
# pixels_mean.to_csv(f'{output_folder}pixel_summary.csv')

# ------------------------visualise------------------------
for_plotting = ratio.copy()
for_plotting['aggregate_cell'] = for_plotting['aggregate_cell'].fillna(0)
for_plotting['sample_key'] = for_plotting['treatment'] +' '+ for_plotting['chaperone']


color_dict = {1: 'rebeccapurple', 0: 'darkorange'}

fig, ax = plt.subplots()
sns.swarmplot(x='sample_key', y='nuc-cyto_ratio', data=for_plotting, hue='aggregate_cell', palette=color_dict, dodge=True)
plt.legend(title='Aggregate cell')
plt.xlabel('Sample')
plt.ylabel('Mean intensity ratio (Nucleus/Cytoplasm)')
plt.show()
