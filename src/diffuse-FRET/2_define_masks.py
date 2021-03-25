
import napari
import numpy as np
import matplotlib.pyplot as plt
import os, shutil
import collections
import skimage.io
from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label
from skimage.morphology import closing, square, remove_small_objects

from loguru import logger

image_folder = f'results/example_diffuse-FRET/initial_cleanup/'
mask_folder = f'results/example_diffuse-FRET/cellpose_masking/'
output_folder = f'results/example_diffuse-FRET/napari_masking/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def filter_masks(image_stack, image_name, mask_stack):

    barnase = mask_stack[0, :, :].copy()
    # Setting nuc masks and agg masks to single value so they are same colour
    # easier to identify when editing masks
    # user will relabel with matching numbers to corresponding cells
    nuc_mask = np.where(mask_stack[1, :, :].copy() != 0, 200, 0)
    incl = np.where(mask_stack[2, :, :].copy() != 0, 100, 0)

    with napari.gui_qt():
        # create the viewer and add the image
        viewer = napari.view_image(image_stack, name='image_stack')
        # add the labels
        viewer.add_labels(barnase, name='barnase')
        viewer.add_labels(incl, name='aggregates')
        viewer.add_labels(nuc_mask, name='mask_features')
        
        """
        - Select the cell layer and using the fill tool set to 0, remove all unwanted cells.
        - Repeat with the inclusions and mask_features layer.
        - Next, select the cell layer and reassign each cell of interest sequentially using the fill tool so that cells are numbered 1 --> n.
        - Repeat with inclusions and mask_features layer, such that inclusion and feature labels correspond to the cell number of interest.
        - Finally, using the brush tool add or adjust any additional features (e.g. Barnase inclusions not associated with incl should be added to the mask_features layer to be removed from diffuse Barnase).
        """
    # collect shapes from inclusions into labels --> can this be coloured easily?

    np.save(f'{output_folder}{image_name}_mask.npy', np.stack([barnase, incl, nuc_mask]))
    logger.info(f'Processed {image_name}. Mask saved to {output_folder}{image_name}')

    return np.stack([barnase, incl, nuc_mask])


# --------------Initialise file list--------------

# reading in all images, and transposing to correct dimension of array
file_list = [filename for filename in os.listdir(image_folder) if '.tif' in filename]
images = {filename.replace('.tif', ''): skimage.io.imread(f'{image_folder}{filename}') for filename in file_list}

# with napari.gui_qt():
#     viewer = napari.view_image(images.values()[0])

# ----------read in masks----------
wholecell = np.load(f'{mask_folder}cellpose_masks.npy')
nucleus = np.load(f'{mask_folder}cellpose_nuclei.npy')
incl = np.load(f'{mask_folder}cellpose_inclusions.npy')

# interleave masks to create single stack per image
raw_masks = {}
for x, image_name in (enumerate(images.keys())):
    raw_masks[image_name] = np.stack([wholecell[x, :, :], nucleus[x, :, :], incl[x, :, :]])

# Manually filter masks, label according to grouped features (i.e. one cell, nucleus (optional) and inclusion per cell of interest, with individual labels)
filtered_masks = {}
for image_name, image_stack in images.items():
    mask_stack = raw_masks[image_name].copy()
    filtered_masks[image_name] = filter_masks(image_stack, image_name, mask_stack)

# # to reprocess individual images:
# filtered_masks = {}
# images_to_process = ['WT_compiled_8']
# for image_name, image_stack in images.items():
#     if image_name in images_to_process:
#         mask_stack = raw_masks[image_name].copy()
#         # remove existing masks
#         if os.path.exists(f'{output_folder}{image_name}/'):
#             shutil.rmtree(f'{output_folder}{image_name}/') 
#         filtered_masks[image_name] = filter_masks(image_stack, image_name, mask_stack)

# # To reload previous masks for per-cell extraction
filtered_masks = {masks.replace('_mask.npy', ''): np.load(f'{output_folder}{masks}') for masks in os.listdir(f'{output_folder}') if '.npy' in masks}

# For each set of masks, separate according to cell number
final_masks = {}
for image_name, image in images.items():
    image_name
    mask_stack = filtered_masks[image_name].copy()
    # plt.imshow(cyto_mask+nuc_mask)
    for cell_number in np.unique(mask_stack[0, :, :]):
        logger.info(cell_number)
        if cell_number > 0:
            whole_cell = np.where(mask_stack[0, :, :] == cell_number, 1, 0)
            mask = np.where(mask_stack[2, :, :] == cell_number, 1, 0)
            aggregates = np.where(mask_stack[1, :, :] == cell_number, 1, 0)
            barnase = np.where(mask == 0, whole_cell, 0)
            barnase = np.where(aggregates == 0, barnase, 0)
            final_masks[(image_name, cell_number)] = np.stack([barnase, aggregates, np.where(mask == 0, whole_cell, 0)])

# ------------------save arrays------------------
for (image_name, cell_number), array_stack in final_masks.items():

    #create folder for each image output
    if not os.path.exists(f'{output_folder}{image_name}/'):
        os.makedirs(f'{output_folder}{image_name}/')
    
    # save associated cell mask arrays
    np.save(f'{output_folder}{image_name}/cell_{int(cell_number)}.npy', array_stack)