import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io
import functools

from loguru import logger

from utilities.file_handling import df_to_excel

logger.info('Import OK')

# define location parameters
input_folder = f'results/chaperone_localisation/pixel_collection/'
output_folder = f'results/chaperone_localisation/summary_calculations/'


if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# read in calculated pixel data
file_list = [filename for filename in os.listdir(input_folder) if '.csv' in filename]
pixels = {filename.replace('.csv', ''): pd.read_csv(f'{input_folder}{filename}') for filename in file_list}

# generate summary df for mask, channel of interest
pixels_compiled = pd.concat(pixels.values())
pixels_compiled.drop([col for col in pixels_compiled.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# determine mean per object (cell, cytoplasm, nuclei) for each channel for each sample
quant_cols = [col for col in pixels_compiled.columns.tolist() if 'channel_' in col]
summary = pixels_compiled.groupby(['image_name', 'cell_number', 'mask_type']).agg({col: ['mean', 'std'] for col in quant_cols})
# fix column names to remove multi-index
summary.columns = ['_'.join(col).strip() for col in summary.columns.values]
summary.reset_index(inplace=True)

# Assign per-treatment identifiers
summary[['chaperone', 'mutant', 'replicate']] = summary['image_name'].str.split('_', expand=True)

# generate nucleus vs cytoplasm ratio
quant_cols = [col for col in summary.columns.tolist() if 'mean' in col]
info_cols = ['image_name', 'cell_number', 'chaperone', 'mutant', 'replicate']
cytoplasm = summary[summary['mask_type'] == 'cytoplasm'].copy().set_index(info_cols)[quant_cols]
nucleus = summary[summary['mask_type'] == 'nucleus'].copy().set_index(info_cols)[quant_cols]

ratio = cytoplasm / nucleus


# save to excel
df_to_excel(
    output_path=f'{output_folder}summary_calculations.xlsx', 
    sheetnames=['summary', 'cyto-nuc_ratio'], 
    data_frames=[summary, ratio.reset_index()])

# visualise example channel 2
info_cols = ['image_name', 'cell_number', 'chaperone', 'mutant', 'replicate']
for_plotting = ratio.reset_index().copy()

color_dict = {'Hsp40': 'firebrick', 'Hsp70': 'darkorange'}

fig, ax = plt.subplots()
sns.boxplot(x='mutant', y=f'channel_2_mean', data=for_plotting.groupby(['chaperone', 'mutant']).mean().reset_index(), hue='chaperone')
sns.stripplot(x='mutant', y=f'channel_2_mean', data=for_plotting, hue='chaperone', palette=color_dict, alpha=0.7, dodge=True)
# Get the handles and labels.
handles, labels = ax.get_legend_handles_labels()
# When creating the legend, only use the second half of the elements to use those from scatter
l = plt.legend(handles[int(len(handles)/2):], labels[int(len(handles)/2):], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.ylabel('Mean ratio (cytoplasm / nucleus)')
plt.xlabel('Variant')
plt.show()