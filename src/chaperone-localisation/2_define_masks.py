
import napari
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
import collections
import skimage.io
from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label
from skimage.morphology import closing, square, remove_small_objects
from loguru import logger


image_folder = f'results/chaperone_localisation/initial_cleanup/'
mask_folder = f'results/chaperone_localisation/cellpose/'
output_folder = f'results/chaperone_localisation/napari_masking/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def filter_masks(image_stack, image_name, mask_stack):

    barnase = mask_stack[0, :, :].copy()
    nuc_mask = mask_stack[1, :, :].copy()
    htt_inc = np.where(mask_stack[2, :, :].copy() != 0, 100, 0)

    with napari.gui_qt():
        # create the viewer and add the image
        viewer = napari.view_image(image_stack, name='image_stack')
        # add the labels
        viewer.add_labels(barnase, name='barnase')
        viewer.add_labels(htt_inc, name='aggregates')
        viewer.add_labels(nuc_mask, name='mask_features')

        """
        - Select the cell layer and using the fill tool set to 0, remove all unwanted cells.
        - Repeat with the inclusions and mask_features layer.
        - Next, select the cell layer and reassign each cell of interest sequentially using the fill tool so that cells are numbered 1 --> n.
        - Repeat with inclusions and mask_features layer, such that inclusion and feature labels correspond to the cell number of interest.
        - Finally, using the brush tool add or adjust any additional features (e.g. Barnase inclusions not associated with Htt should be added to the mask_features layer to be removed from diffuse Barnase).
        """
    # collect shapes from inclusions into labels --> can this be coloured easily?

    np.save(f'{output_folder}{image_name}_mask.npy',
            np.stack([barnase, htt_inc, nuc_mask]))
    logger.info(
        f'Processed {image_name}. Mask saved to {output_folder}{image_name}')

    return np.stack([barnase, htt_inc, nuc_mask])


# --------------Initialise file list--------------

# reading in all images, and transposing to correct dimension of array
file_list = [filename for filename in os.listdir(
    image_folder) if '.tif' in filename]
images = {filename.replace('.tif', ''): skimage.io.imread(
    f'{image_folder}{filename}').transpose(2, 0, 1) for filename in file_list}

with napari.gui_qt():
    viewer = napari.view_image(list(images.values())[0][:, :, :])

# ----------read in masks----------
wholecell = np.load(f'{mask_folder}cellpose_masks.npy')
nucleus = np.load(f'{mask_folder}cellpose_nuclei.npy')
htt = np.load(f'{mask_folder}cellpose_inclusions.npy')

# interleave masks to create single stack per image
raw_masks = {}
for x, image_name in (enumerate(images.keys())):
    raw_masks[image_name] = np.stack(
        [wholecell[x, :, :], nucleus[x, :, :], htt[x, :, :]])

# Manually filter masks, label according to grouped features (i.e. one cell, nucleus (optional) and inclusion per cell of interest, with individual labels)
filtered_masks = {}
for image_name, image_stack in images.items():
    mask_stack = raw_masks[image_name].copy()
    filtered_masks[image_name] = filter_masks(
        image_stack, image_name, mask_stack)


# --------------------- To reload previous masks for per-cell extraction---------------------
filtered_masks = {masks.replace('_mask.npy', ''): np.load(
    f'{output_folder}{masks}') for masks in os.listdir(f'{output_folder}') if '.npy' in masks}

# For each set of masks, separate according to cell number
final_masks = {}
for image_name, image in images.items():
    image_name
    mask_stack = filtered_masks[image_name].copy()
    # plt.imshow(cyto_mask+nuc_mask)
    for cell_number in np.unique(mask_stack[0, :, :]):
        logger.info(cell_number)
        if cell_number > 0: # background is currently 0, cells are numbered sequentially from 1 -> n
            # select individual cell where the mask is equal to that cell number, replace that cell number with 1's and fill the rest of the mask with 0
            whole_cell = np.where(mask_stack[0, :, :] == cell_number, 1, 0)

            # where the whole cell mask is equal to 1, get the nucleus pixels from nuc_mask and fill the rest of the mask with 0
            nucleus = np.where(whole_cell == 1, mask_stack[2, :, :], 0)
            # where the nucleus mask is anything other than 0, change it to 1, then fill the rest of the mask with 0
            nucleus = np.where(nucleus != 0, 1, 0)

            # repeat steps as above, but for the agg mask
            aggregates = np.where(whole_cell == 1, mask_stack[1, :, :], 0)
            aggregates = np.where(aggregates != 0, 1, 0)


            # where the nucleus mask is 0, get the whole cell mask (remembering that the wc mask is 1 where cell is, 0 where it's not), then fill the rest of the mask with 0
            cytoplasm = np.where(nucleus == 0, whole_cell, 0)
            plt.imshow(cytoplasm)

            
            cytoplasm = np.where(aggregates == 0, cytoplasm, 0) # exclude aggregate from cyto
            nucleus = np.where(aggregates == 0, nucleus, 0) # exclude aggregate from cyto


            final_masks[(image_name, cell_number)] = np.stack(
                [cytoplasm, aggregates, nucleus])



# ------------------save arrays------------------
for (image_name, cell_number), array_stack in final_masks.items():

    #create folder for each image output
    if not os.path.exists(f'{output_folder}{image_name}/'):
        os.makedirs(f'{output_folder}{image_name}/')

    # save associated cell mask arrays
    np.save(f'{output_folder}{image_name}/cell_{int(cell_number)}.npy', array_stack)