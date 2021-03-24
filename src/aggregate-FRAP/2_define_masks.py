
import napari
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import collections
import skimage.io
from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label
from skimage.morphology import closing, square, remove_small_objects

from loguru import logger

input_folder = f'results/aggregate-FRAP/initial_cleanup/'
output_folder = f'results/aggregate-FRAP/napari_masking/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def mask_per_stack(image_stack, image_name, num_pre=5, num_bleach=30):
    """Define bleaching ROI via thresholding, enable editable non-bleached and background ROIs of equal size. ROI positions are then mapped across all timepoints.

    Parameters
    ----------
    image_stack : array
        original image stack where z-dimension are timepoints
    image_name : str
        name of image to be processed
    num_pre : int, optional
        number of frames taken pre-bleach, by default 5
    num_bleach : int, optional
        number of frames used for bleaching, by default 30

    Returns
    -------
    tuple(df, dict)
        coords: df mapping x, y pixels inside the thresholded bleach ROI
        timepoints: dict mapping each timepoint to ROI mask array containing background, non-bleached and bleached masks
    """
    # Threshold fist bleach image to create bleach mask
    bleach_image = image_stack[:, :, num_pre+1]
    # plt.imshow(bleach_image)
    # apply threshold
    thresh = threshold_otsu(bleach_image)
    bw = closing(bleach_image > thresh, square(4))
    # label image regions
    bleach_ROIs = label(bw)

    # calculate radius and centroid for circle shape
    # get list of positions for each ROI
    pixel_array = np.where(bleach_ROIs != 0, bleach_ROIs, np.nan)
    # plt.imshow(pixel_array)
    coords = pd.DataFrame(pixel_array).unstack().reset_index().dropna()
    coords.columns = ['y', 'x', 'label']
    # calculate position info for shape data 
    # centroid is mean x+y positions
    centroids = coords.groupby('label').mean()
    # radius is max - min /2
    centroids['radii'] = ((coords.groupby('label').max() - coords.groupby('label').min()) / 2).min(axis=1)

    masks = {}

    for roi_label, roi in centroids.iterrows():
        bleach_roi_shape = np.array([[roi['x'], roi['y']], [roi['radii'], roi['radii']]])
        background_roi_shape = np.array([[roi['radii']+1, roi['radii']+1], [roi['radii'], roi['radii']]])
        nonbleach_roi_shape = np.array([[roi['x'], roi['y']], [roi['radii'], roi['radii']]])

        with napari.gui_qt():
            # create the viewer and add the image
            viewer = napari.view_image(image_stack.transpose(2, 0, 1), name='image_stack')
            # add the labels
            viewer.add_labels(bleach_ROIs, name='segmentation')
            # add shapes
            b_layer = viewer.add_shapes(bleach_roi_shape, shape_type='ellipse', edge_width=1, name='bleach')
            nb_layer = viewer.add_shapes(nonbleach_roi_shape, shape_type='ellipse', edge_width=1, name='non-bleach')
            bg_layer = viewer.add_shapes(background_roi_shape, shape_type='ellipse', edge_width=1, name='background')

        bg_mask = bg_layer.to_labels((256, 256))
        nb_mask = nb_layer.to_labels((256, 256))
        b_mask = b_layer.to_labels((256, 256))

        timepoints = {
            timepoint: np.array([bg_mask, nb_mask, b_mask])
            for timepoint in range(image_stack.shape[2])
        }

        masks[f'{image_name}_{int(roi_label)}'] = timepoints

    return coords, masks


def mask_per_timepoint(image_stack, image_name, num_pre=5, num_bleach=30, visualise_timepoints=False):
    """Spawn editable ROIs for individual timepoints, which are then carried over subsequent timepoints until next edit point. Useful for cells that are moving, or with long per-frame intervals where the position of the thresholded bleaching ROI is insufficient to set the bleached ROI over the recovery period.

    Parameters
    ----------
    image_stack : array
        original image stack where z-dimension are timepoints
    image_name : str
        name of image to be processed
    num_pre : int, optional
        number of frames taken pre-bleach, by default 5
    num_bleach : int, optional
        number of frames used for bleaching, by default 30
    visualise_timepoints : OrderedDict, optional
        Dictionary mapping number of frames post-bleach to how often images should be displayed for ROI editing
        e.g. visualise_timepoints={30: 5, 20: 2, 10: 1} visualises every 5 frames from 35-65,  every 5 frames from 65-85 and every frame from 85-95 assuming default num_pre and num_bleach. Any unaccounted for frames from end of visualise_timepoints to total frames will be assigned identical ROI to the last visualise_timepoints frame, so when specifying this method it is recommended to account for all frames.
        By default False to visualise all timepoints.

    Returns
    -------
    tuple(df, dict)
        coords: df mapping x, y pixels inside the thresholded bleach ROI
        timepoints: dict mapping each timepoint to ROI mask array containing background, non-bleached and bleached masks
    """
    # Threshold fist bleach image to create bleach mask
    bleach_image = image_stack[:, :, num_pre+1]
    # plt.imshow(bleach)
    # apply threshold
    thresh = threshold_otsu(bleach_image)
    bw = closing(bleach_image > thresh, square(4))
    # label image regions
    bleach_ROIs = label(bw)

    # calculate radius and centroid for circle shape
    # get list of positions for each ROI
    pixel_array = np.where(bleach_ROIs != 0, bleach_ROIs, np.nan)
    # plt.imshow(pixel_array)
    coords = pd.DataFrame(pixel_array).unstack().reset_index().dropna()
    coords.columns = ['y', 'x', 'label']
    # calculate position info for shape data 
    # centroid is mean x+y positions
    centroids = coords.groupby('label').mean()
    # radius is max - min /2
    centroids['radii'] = ((coords.groupby('label').max() - coords.groupby('label').min()) / 2).min(axis=1)

    if visualise_timepoints:
        # Visualise timepoints should be provided in as dict(number of frames, steps)
        frames = [(num_pre+num_bleach)] + list(visualise_timepoints.keys())
        frames = np.cumsum(frames)
        frames_to_view = [list(np.arange(frames[x], frames[x + 1] + 1, step))for x, step in enumerate(visualise_timepoints.values())        ]

        frames_to_view = list({item for sublist in frames_to_view for item in sublist})

    masks = {}

    for roi_label, roi in centroids.iterrows():
        bleach_roi = np.array([[roi['x'], roi['y']], [roi['radii'], roi['radii']]])
        background_roi = np.array([[roi['radii']+1, roi['radii']+1], [roi['radii'], roi['radii']]])
        nonbleach_roi = np.array([[roi['x'], roi['y']], [roi['radii'], roi['radii']]])

        # apply ROIs to all post-bleached images
        timepoints = {}
        for timepoint in range((num_pre+num_bleach), image_stack.shape[2]):
            if timepoint in frames_to_view:
                with napari.gui_qt():
                    # create the viewer and add the image
                    viewer = napari.view_image(image_stack[:, :, timepoint], name='image')
                    # add the labels
                    viewer.add_labels(bleach_ROIs, name='segmentation')
                    # add shapes
                    b_layer = viewer.add_shapes(bleach_roi, shape_type='ellipse', edge_width=1, name=f'bleach_{timepoint}')
                    nb_layer = viewer.add_shapes(nonbleach_roi, shape_type='ellipse', edge_width=1, name=f'non-bleach_{timepoint}')
                    bg_layer = viewer.add_shapes(background_roi, shape_type='ellipse', edge_width=1, name=f'background_{timepoint}')

            bg_mask = bg_layer.to_labels((256, 256))
            nb_mask = nb_layer.to_labels((256, 256))
            b_mask = b_layer.to_labels((256, 256))
            timepoints[timepoint] = np.array([bg_mask, nb_mask, b_mask])

            bleach_roi = b_layer.data
            nonbleach_roi = nb_layer.data
            background_roi = bg_layer.data

        # apply first post-bleach ROIs to all pre-bleach and bleach images
        for timepoint in range((num_pre+num_bleach)):
            timepoints[timepoint] = timepoints[(num_pre+num_bleach)]  

        masks[f'{image_name}_{int(roi_label)}'] = timepoints

    return coords, masks


# --------------Initialise file list--------------
# read in pre-stacked arrays
file_list = [filename for filename in os.listdir(input_folder) if '.npy' in filename]

# reading in all images, and transposing to correct dimension of array
images = {filename.replace('.npy', ''): np.load(f'{input_folder}{filename}') for filename in file_list}

# with napari.gui_qt():
#    viewer = napari.view_image(images['example_1'].transpose(2, 0, 1))

# ---------to process a subset of images---------
images_to_process = []

if len(images_to_process) < 1:
    images_to_process = images.keys()
    logger.info('Processing all images')

# ----------generate masks for each ROI----------
for image_name in images_to_process:
    image = images[image_name]
    # coords, mask = mask_per_stack(image_stack=image, image_name=image_name, num_pre=5, num_bleach=30)

    visualise_timepoints={
    30: 5, # for 30 frames post-beach, show me every 5 frames
    20: 2, # for the next 20 frames, show me every 2 frames
    10: 1, # for the next 10 frames, show me every frame
    }
    coords, mask = mask_per_timepoint(image_stack=image, image_name=image_name, num_pre=5, num_bleach=30, visualise_timepoints=visualise_timepoints)

    # -----------------save arrays-----------------
    for roi_name, timepoints in mask.items():
        if not os.path.exists(f'{output_folder}{roi_name}/'):
            os.makedirs(f'{output_folder}{roi_name}/')

        for timepoint, array_stack in timepoints.items():
            # save associated arrays
            np.save(f'{output_folder}{roi_name}/{timepoint}.npy', array_stack)
