import os
import numpy as np
import matplotlib.pyplot as plt
import skimage.io
from cellpose import models
from cellpose import plot
import collections


input_folder = f'results/chaperone_localisation/initial_cleanup/'
output_folder = f'results/chaperone_localisation/cellpose_masking/'

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
img_names = [filename for filename in file_list]

# -----------------------Complete cellpose with cytoplasm channel---------------------------------
#Channel 0: Hoerscht
#Channel 1: BF
#Channel 2: Antibody - this looks ok, and certainly better than brightfield. Need to be careful for those that are clumped together - may need to consider Bf

# collecting only channel 2's for cytoplasm
cytoplasm_images = [image[:, :, 2] for image in imgs]
plt.imshow(cytoplasm_images[0])

# Apply cellpose then visualise
masks, flows, styles, diams = apply_cellpose(cytoplasm_images, image_type='cyto', diameter=300)
visualise_cell_pose(cytoplasm_images, masks, flows, channels=[0, 0])

# -----------------------Complete cellpose with nuclei channel---------------------------------
# collecting only channel 0's for cytoplasm
nuclei_images = [image[:, :, 0] for image in imgs]
# plt.imshow(cytoplasm_images[0]+nuclei_images[0])

# Apply cellpose then visualise
nuc_masks, nuc_flows, nuc_styles, nuc_diams = apply_cellpose(nuclei_images, diameter=200, image_type='nuclei')
visualise_cell_pose(nuclei_images, nuc_masks, nuc_flows, channels=[0, 0])

# save associated cell mask arrays
np.save(f'{output_folder}cellpose_masks.npy', masks)
np.save(f'{output_folder}cellpose_nuclei.npy', nuc_masks)

# optional - edit nuclei masks?
# nuc_masks = segmentation_editor(images, masks)

