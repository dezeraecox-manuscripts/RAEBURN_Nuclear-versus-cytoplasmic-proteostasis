import os
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
import skimage.io
import functools


from GEN_Utils.FileHandling import df_to_excel
from loguru import logger

logger.info('Import OK')

# define location parameters
image_folder = f'results/chaperone_localisation/initial_cleanup/'
mask_folder = f'results/chaperone_localisation/napari_masking/'
output_folder = f'results/chaperone_localisation/pixel_collection/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# define function to collect pixel location and intensity for a given mask, image combination
def pixel_collector(image_array, mask, mask_type=None, visualise=False, size=(1024, 1024)):
    pixel_array = np.where(mask != 0, image_array, np.nan)
    # plt.imshow(pixel_array)
    coords = pd.DataFrame(pixel_array).unstack().reset_index().dropna()
    coords.columns = ['y', 'x', 'intensity']

    if mask_type != None:
        coords['mask_type'] = mask_type

    if visualise:
        # test visualisation, compare to plt.show
        fig, ax = plt.subplots(figsize=(20, 20))
        sns.scatterplot(coords['x'], coords['y'], hue=coords['intensity'], palette='magma_r', size=0.5,linewidth=0, alpha = 0.7)
        plt.ylim(size[0], 0)
        plt.xlim(0, size[1])

    return coords

# --------------Initialise file lists--------------
# reading in all images, and transposing to correct dimension of array
images = { image_name.replace('.tif', ''): skimage.io.imread(f'{image_folder}{image_name}') for image_name in os.listdir(f'{image_folder}') if '.tif' in image_name}

# read in masks - remember that stack format is [barnase, aggregates]
# fails try/except if no masks found therefore skip that image
masks = {}
for image_name in images.keys():
    logger.info(f'Processing {image_name}')
    try:
        masks[f'{image_name}'] = {cell.replace('.npy', ''): np.load(f'{mask_folder}{image_name}/{cell}') for cell in os.listdir(f'{mask_folder}{image_name}/')}
        logger.info(f'Masks loaded for {len(masks[f"{image_name}"].keys())} cells')
    except:
        logger.info(f'{image_name} not processed as no mask found')


# Example napari visualisation to test matching mask array to image
import napari
image_test_name = '72Q Httex1_HSP40_Series003'
cyto_mask = np.zeros((1024, 1024))
nuc_mask = np.zeros((1024, 1024))
agg_mask = np.zeros((1024, 1024))
with napari.gui_qt():
    viewer = napari.Viewer()
    viewer.add_image(images[image_test_name].transpose(2, 0, 1), name='raw_image')
    for cell in masks[f'{image_test_name}'].keys():
        cyto_mask = np.where(masks[f'{image_test_name}'][cell][0, :, :] == 1, int(cell.split('_')[-1]), cyto_mask)
        agg_mask = np.where(masks[f'{image_test_name}'][cell][1, :, :] == 1, int(cell.split('_')[-1]), agg_mask)
        nuc_mask = np.where(masks[f'{image_test_name}'][cell][2, :, :] == 1, int(cell.split('_')[-1]), nuc_mask)

    viewer.add_labels(cyto_mask, name=f'cytoplasm')
    viewer.add_labels(agg_mask, name=f'aggregate')
    viewer.add_labels(nuc_mask, name=f'nucleus')



# ---------------collect pixel information---------------
pixel_information = {}
for image_name in masks.keys():
    logger.info(f'Processing {image_name}')
    image = images[image_name].copy()
    pixels = []
    for channel in range(image.shape[2]):
        image_array = image[:, :, channel].copy()
        for cell, mask_stack in masks[image_name].items():
            # for each cell, collect barnase and aggregate pixels
            cell_pixels = pd.concat([pixel_collector(image_array, mask_stack[i, :, :], visualise=False, mask_type=mask_type) for i, mask_type in enumerate(['cytoplasm', 'aggregate', 'nucleus'])])
            # add identifiers
            cell_pixels['channel'] = channel
            cell_pixels['cell'] = f'{image_name}_{cell}'
            pixels.append(cell_pixels)
    pixels = pd.concat(pixels)
    pixel_information[image_name] = pixels
logger.info('Completed pixel collection')

# save to excel (although this will likely be unopenable if more than a few images) and csv
# df_to_excel(output_path=f'{output_folder}pixel_information.xlsx', sheetnames=list(pixel_information.keys()), data_frames=list(pixel_information.values()))
saved = [df.to_csv(f'{output_folder}{image_name}.csv') for image_name, df in pixel_information.items()]
