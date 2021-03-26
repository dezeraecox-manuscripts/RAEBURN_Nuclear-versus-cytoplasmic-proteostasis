import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io
import functools

from utilities.pixel_operations import pixel_collector

from loguru import logger

logger.info('Import OK')

# define location parameters
image_folder = f'results/aggregate-FRAP/initial_cleanup/'
mask_folder = f'results/aggregate-FRAP/napari_masking/'
output_folder = f'results/aggregate-FRAP/pixel_collection/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# --------------Initialise file lists--------------
mask_list = [folder for folder in os.listdir(mask_folder) if '.DS' not in folder]
image_names = list({'_'.join(folder.split('_')[:-1]) for folder in mask_list})

# reading in all images, and transposing to correct dimension of array
images = {filename.replace('.npy', ''): np.load(f'{image_folder}{filename}.npy') for filename in image_names}


# read in masks, collect timepoints
# - remember that stack format is [background, nonbleach, bleach]
masks = {}
for roi_name in mask_list:
    logger.info(f'Processing {roi_name}')
    try:
        masks[roi_name] = {timepoint.replace('.npy', ''): np.load(f'{mask_folder}{roi_name}/{timepoint}') for timepoint in os.listdir(f'{mask_folder}{roi_name}/')}
        logger.info(f'Masks loaded for {len(masks[roi_name].keys())} timepoints')
    except:
        logger.info(f'{image_name} not processed as no mask found')


# # Example napari visualisation to test matching mask array to image
# import napari
# image_test_name = '4Y_3'
# with napari.gui_qt():
#     viewer = napari.Viewer()
#     viewer.add_image(images[image_test_name][:, :, 0], name='raw_image')
#     viewer.add_labels(masks[f'{image_test_name}_1']['0'][0, :, :], name='background')
#     viewer.add_labels(masks[f'{image_test_name}_1']['0'][1, :, :], name='nonbleach')
#     viewer.add_labels(masks[f'{image_test_name}_1']['0'][2, :, :], name='bleach')

# ---------------collect pixel information---------------
pixel_information = {}
for roi_name in mask_list:
    logger.info(f'Processing {roi_name}')
    image_name = '_'.join(roi_name.split('_')[:-1])
    image = images[image_name]
    timepoints = []
    for timepoint in range(image.shape[2]):
        # logger.info(f'Processing timepoint {timepoint} pixels')
        # collect only one timepoint
        image_array = image[:, :, timepoint]
        mask_array = masks[roi_name][f'{timepoint}']
        # for each timepoint, collect background, bleach and non-bleach pixels
        roi_pixels = pd.concat([pixel_collector(image_array, mask_array[i, :, :], visualise=False, mask_type=mask_type) for i, mask_type in enumerate(['background', 'nonbleach', 'bleach'])])
        # add identifiers
        roi_pixels['timepoint'] = timepoint
        timepoints.append(roi_pixels)
    timepoints = pd.concat(timepoints)
    timepoints['roi_name'] = roi_name

    pixel_information[roi_name] = timepoints
logger.info('Completed pixel collection')

# save to csv
saved = [df.to_csv(f'{output_folder}{roi_name}.csv') for roi_name, df in pixel_information.items()]
