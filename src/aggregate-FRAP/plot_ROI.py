import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger

logger.info('Import OK')

input_path = 'results/aggregate-FRAP/summary_calculations/FRAP_summary.xlsx'
image_folder = 'results/aggregate-FRAP/initial_cleanup/'
output_folder = 'results/aggregate-FRAP/plot_ROI/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# read in pixel information
pixel_summary = pd.read_excel(f'{input_path}', sheet_name='compiled')
pixel_summary.drop([col for col in pixel_summary.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Calculate centre of each ROI
pixel_summary = pixel_summary[pixel_summary['timepoint'] == 0]
roi_centroid = pixel_summary.groupby(['roi_name', 'mask_type']).mean().reset_index()
roi_centroid['image_name'] = roi_centroid['roi_name'].str.split('_').str[:-1].str.join('_')
pixel_summary['image_name'] = pixel_summary['roi_name'].str.split('_').str[:-1].str.join('_')

# read in images
file_list = [filename for filename in os.listdir(image_folder) if '.npy' in filename]

# reading in all images, and collecting t0
images = {filename.replace('.npy', ''): np.load(f'{image_folder}{filename}') for filename in file_list}

# Generate plots for each image_name
for image_name, df in roi_centroid.groupby('image_name'):
    fig, ax = plt.subplots()
    ax.imshow(images[image_name][:, :, 35], cmap='Greys_r')
    ax.scatter(data=df, x='x', y='y', color='red')
    for roi_name, roi_type, x, y in df[['roi_name', 'mask_type', 'x', 'y']].values:
        ax.annotate(f'{roi_name}_{roi_type}', (x+10, y+10), c='red', fontsize=6)
    plt.savefig(f'{output_folder}{image_name}.png')