from skimage.morphology import closing, square, remove_small_objects
from skimage.measure import label
from skimage.segmentation import clear_border
from skimage.filters import threshold_otsu
import os
import numpy as np
import matplotlib.pyplot as plt
import skimage.io
from cellpose import models
from cellpose import plot
import collections


input_folder = f'results/example_diffuse-FRET/initial_cleanup/'
output_folder = f'results/example_diffuse-FRET/cellpose_masking/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

def apply_cellpose(images, image_type='cyto', channels=[0,0], diameter=None):
    """Apply model to list of images. Returns masks, flows, styles, diams.
    - model type is 'cyto' or 'nuclei'
    - define CHANNELS to run segementation on (grayscale=0, R=1, G=2, B=3) where channels = [cytoplasm, nucleus]. If NUCLEUS channel does not exist, set the second channel to 0
    """
    model = models.Cellpose(model_type=image_type)
    masks, flows, styles, diams = model.eval(images, diameter=diameter, channels=channels)
    return masks, flows, styles, diams

def visualise_cell_pose(images, masks, flows, channels=[0,0]):
    """Display cellpose results for each image
    """
    for image_number, image in enumerate(images):
        maski = masks[image_number]
        flowi = flows[image_number][0]

        fig = plt.figure(figsize=(12,5))
        plot.show_segmentation(fig, image, maski, flowi, channels=channels)
        plt.tight_layout()
        plt.show()

def edge_filter(mask):
    """Collect boundary pixel values for all edges, return unique values
    which correspond to cells that are touching/over the edge boundaries""" 
    size = mask.shape[0]
    edge_1_cells = set(list(mask[0:1, 0:size].flatten()))
    edge_2_cells = set(list(mask[0:size, (size-1):size].flatten()))
    edge_3_cells = set(list(mask[(size-1):size, 0:size].flatten()))
    edge_4_cells = set(list(mask[0:size, 0:1].flatten()))
    return edge_1_cells | edge_2_cells | edge_3_cells | edge_4_cells

def size_filter(mask, lower_size=1500, upper_size=10000):
    """Collect cells that are outside the cell size bounds as those to
    be excluded""" 
    cell_size = dict(collections.Counter(mask.flatten()))
    bs_cells = [
        cell_number
        for cell_number, cell_size in cell_size.items()
        if cell_size < lower_size or cell_size > upper_size
    ]

    return set(bs_cells)

# --------------------------------------Initialise file list--------------------------------------

file_list = [filename for filename in os.listdir(input_folder) if '.tif' in filename]

# reading in all channels for each image, and transposing to correct dimension of array
imgs = [skimage.io.imread(f'{input_folder}{filename}').transpose(1, 2, 0) for filename in file_list]

# clean filenames
img_names = [filename.replace('.tif', '') for filename in file_list]

# -----------------------Complete cellpose with cytoplasm channel---------------------------------
# channel 0: Venus ----> use this for making masks
# channel 1: Bright field
# channel 2: mTFP
# channel 3: FRET
# channel 4: inclusions

# collecting only channel 0's for masking
cytoplasm_images = [image[:, :, 0] for image in imgs]
plt.imshow(cytoplasm_images[0])

# Apply cellpose then visualise
masks, flows, styles, diams = apply_cellpose(cytoplasm_images, image_type='cyto', diameter=50)
visualise_cell_pose(cytoplasm_images, masks, flows, channels=[0, 0])

# -----------------------If NES image, use inversion of venus channel to define nuclei---------------------------------
nuc_masks, nuc_flows, nuc_styles, nuc_diams = apply_cellpose([65000 - array for array in cytoplasm_images], image_type='nuclei', diameter=20)
visualise_cell_pose([65000 - array for array in cytoplasm_images], nuc_masks, nuc_flows, channels=[0, 0])


# -----------------------outline inclusions---------------------------------
incl_images = [image[:, :, 4] for image in imgs]
plt.imshow(incl_images[0])
incl_masks, incl_flows, incl_styles, incl_diams = apply_cellpose(incl_images, image_type='nuclei', diameter=20)
visualise_cell_pose(incl_images, incl_masks, incl_flows, channels=[0, 0])



# save associated cell mask arrays
np.save(f'{output_folder}cellpose_masks.npy', masks)
np.save(f'{output_folder}cellpose_nuclei.npy', nuc_masks)
np.save(f'{output_folder}cellpose_inclusions.npy', incl_masks)
