import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io
import functools

from loguru import logger

logger.info('Import OK')

# define location parameters, assign number of frames for pre- and bleach portions
input_folder = f'results/aggregate-FRAP/pixel_collection/'
output_folder = f'results/aggregate-FRAP/summary_calculations/'

num_prebleach = 5
num_bleach = 30

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def df_to_excel(output_path, sheetnames, data_frames):
    """Saves list of dataframes to a single excel (xlsx) file.
    Parameters
    ----------
    output_path : str
        Full path to which xlsx file will be saved.
    sheetnames : list of str
        descriptive list of dataframe content, used to label sheets in xlsx file.
    data_frames : list of DataFrames
        DataFrames to be saved. List order must match order of names provided in sheetname.
    Returns
    -------
    None.
    """
    if not output_path.endswith('.xlsx'):
        output_path = output_path+'Results.xlsx'
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    # Convert the dataframe to an XlsxWriter Excel object.
    for x in range(0, len(sheetnames)):
        sheetname = sheetnames[x]
        data_frame = data_frames[x]
        data_frame.to_excel(writer, sheet_name=sheetname)
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()


# -----Process dataset-----

# read in calculated pixel data
file_list = [filename for filename in os.listdir(input_folder) if '.csv' in filename]
pixels = {filename.replace('.csv', ''): pd.read_csv(f'{input_folder}{filename}') for filename in file_list}
pixels.update({key: value.drop([col for col in value.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, value in pixels.items()})

# generate summary df for mask, timepoint of interest
pixels_compiled = pd.concat(pixels.values())

# generate mean values for each ROI for each timepoint
pixels_mean = pixels_compiled.copy().groupby(['roi_name', 'mask_type', 'timepoint']).mean().reset_index()

# assign timepoint identifiers - this will allow taking mean for pre-bleach, and removing bleach timepoints
timepoint_map = {}
timepoint_map.update(dict(zip(np.arange(0, num_prebleach), [-1] * num_prebleach)))
timepoint_map.update(dict(zip(np.arange(num_prebleach, num_prebleach+num_bleach), [np.nan] * num_bleach)))
timepoint_map.update(dict(zip(np.arange(num_prebleach+num_bleach, pixels_compiled['timepoint'].max()+1), np.arange(0, pixels_compiled['timepoint'].max()+1))))
pixels_mean['timepoint_map'] = pixels_mean['timepoint'].map(timepoint_map)

# Remove bleach timepoints
pixels_mean.dropna(subset=['timepoint_map'], inplace=True)

# pivot table to make FRAP calculations more simple - note this also automatically takes the mean of the prebleach timepoints
pixels_corrected = pd.pivot_table(pixels_mean, values='intensity', index=['roi_name', 'timepoint_map'], columns=['mask_type'], aggfunc='mean').reset_index()

# generate roi - background corrected ratio
# (note this is currently using the individual backround intensity at each timepoint, previously had used mean over all timepoints)
pixels_corrected['bleach_corrected'] = pixels_corrected['bleach'] - pixels_corrected['background']
pixels_corrected['nonbleach_corrected'] = pixels_corrected['nonbleach'] - pixels_corrected['background']
# Calculate FRAP value
pixels_corrected['FRAP'] = pixels_corrected['bleach_corrected'] / pixels_corrected['nonbleach_corrected']
# normalise to nonbleach at prebleach
pixels_summary = []
for roi_name, df in pixels_corrected.groupby('roi_name'):
    prebleach_FRAP = df[df['timepoint_map'] == -1.0]['FRAP']
    df['FRAP_corrected'] = df['FRAP'] / prebleach_FRAP.values[0]
    pixels_summary.append(df)
pixels_summary = pd.concat(pixels_summary)

# Assign sample identifiers
pixels_summary['mutant'] = pixels_summary['roi_name'].str.split('_').str[0]

# save to excel
df_to_excel(
    output_path=f'{output_folder}FRAP_summary.xlsx',
    sheetnames=['summary', 'compiled'],
    data_frames=[pixels_summary, pixels_mean])