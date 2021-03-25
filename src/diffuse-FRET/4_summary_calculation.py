import os
import pandas as pd

from loguru import logger

logger.info('Import OK')

# define location parameters, assign number of frames for pre- and bleach portions
input_folder = f'results/example_diffuse-FRET/pixel_collection/'
output_folder = f'results/example_diffuse-FRET/summary_calculations/'

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
pixels_compiled.columns = ['_'.join([str(val) for val in x]) if type(x[1]) == int else x[0] for x in pixels_compiled.columns]

# Add label if aggregate inside unmasked (i.e. same compartment)
aggregate_cells = pixels_compiled[pixels_compiled['mask_type'] == 'aggregate']['cell'].unique()
aggregate_pixels = pixels_compiled[pixels_compiled['cell'].isin(aggregate_cells)]
aggregate_labels = {}
for cell, df in aggregate_pixels.groupby('cell'):
    agg_loc = set(tuple(zip(df[df['mask_type'] == 'aggregate']['x'], df[df['mask_type'] == 'aggregate']['y'])))
    barnase_loc = set(tuple(zip(df[df['mask_type'] == 'unmasked']['x'], df[df['mask_type'] == 'unmasked']['y'])))
    logger.info(len(agg_loc.intersection(barnase_loc)) / len(agg_loc))
    aggregate_labels[cell] = (round(len(agg_loc.intersection(barnase_loc)) / len(agg_loc), 2), 'inside' if len(agg_loc.intersection(barnase_loc)) / len(agg_loc) > overlap_threshold else 'outside')

# generate mean values for each ROI for each timepoint
pixels_mean = pixels_compiled.copy().groupby(['cell', 'mask_type']).mean().reset_index()

# assign identifiers
pixels_mean[['mutant', 'target', 'image_number', 'discard', 'cell_number']] = pixels_mean['cell'].str.split('_', expand=True)
pixels_mean.drop('discard', axis=1, inplace=True)
pixels_mean[['overlap', 'agg_location']] = pd.DataFrame(pixels_mean['cell'].map(aggregate_labels).tolist(), index=pixels_mean.index)  
pixels_mean['agg_location'] = ['None' if entry == None else entry for entry in pixels_mean['agg_location']]

# save to csv
pixels_mean.to_csv(f'{output_folder}pixel_summary.csv')
