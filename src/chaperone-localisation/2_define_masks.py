# ------------------------Napari to adjust masks-----------------------------------------------

# visualise each mask with napari, allowing for editing and recapture result into "manual masks"
# create copy of masks array as it is overwritten by actions in the viewer
import napari
import numpy as np
import matplotlib.pyplot as plt
import os
import collections
import skimage.io

input_folder = f'results/chaperone_localisation/initial_cleanup/'
output_folder = f'results/chaperone_localisation/cellpose_masking/'


def segmentation_editor(images, masks):
    
    manual_masks = np.copy(masks)
    for i in range(len(manual_masks)):
        with napari.gui_qt():
            # create the viewer and add the coins image
            viewer = napari.view_image(cytoplasm_images[i], name='cells')
            # add the labels
            labels = viewer.add_labels(manual_masks[i], name='segmentation')

    return manual_masks


# read in numpy masks
masks = np.load(f'{output_folder}cellpose_masks.npy')
nuc_masks = np.load(f'{output_folder}cellpose_nuclei.npy')

# --------------------------------------Initialise file list--------------------------------------

file_list = []
for filename in os.listdir(input_folder):
    if filename.endswith('.tif'):
        file_list.append(filename)
    else:
        continue 
# file_list = [filename for filename in os.listdir(input_folder)]
#with napari.gui_qt():
 #   viewer = napari.view_image(skimage.io.imread(f'{input_folder}{file_list[0]}'))

# reading in all channels for each image, and transposing to correct dimension of array
imgs = [skimage.io.imread(f'{input_folder}{filename}').transpose(1, 2, 0) for filename in file_list]
# clean filenames
img_names = [filename.replace('.tif', '') for filename in file_list]
# collecting only channel 0's for cytoplasm
cytoplasm_images = [image[:, :, 2] for image in imgs]

# adjust masks with Napari
filtered_masks = segmentation_editor(cytoplasm_images, masks)
filtered_masks = {image_number: np.where(mask!=0, mask, np.nan) for image_number, mask in enumerate(filtered_masks)}


final_masks = {}
for image_number in range(len(imgs)):
    cyto_mask = filtered_masks[image_number].copy()
    nuc_mask = nuc_masks[image_number].copy()
    # plt.imshow(cyto_mask+nuc_mask)
    for cell_number in np.unique(cyto_mask[~np.isnan(cyto_mask)]):
        whole_cell = np.where(cyto_mask == cell_number, 1, 0)
        nucleus = np.where(whole_cell == 1, nuc_mask, 0)
        # TODO: to filter nuclei for all within cell, check count nuc_id before + after filtering
        # TODO: to filter nuclei for only one per cell, check in len(unique_nuc) == 1
        cytoplasm = np.where(nucleus == 0, whole_cell, 0)
        nucleus = np.where(nucleus != 0, 1, 0)
        final_masks[(image_number, cell_number)] = np.dstack([whole_cell, nucleus, cytoplasm])

# -----------------------------------------save arrays-----------------------------------------
for (image_number, cell_number), array_stack in final_masks.items():

    #create folder for each image output
    if not os.path.exists(f'{output_folder}{img_names[image_number]}/'):
        os.mkdir(f'{output_folder}{img_names[image_number]}/')
    
    # save associated cell mask arrays
    np.save(f'{output_folder}{img_names[image_number]}/cell_{int(cell_number)}.npy', array_stack)
