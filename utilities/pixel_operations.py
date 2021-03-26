import os
import numpy as np
import pandas as pd
import skimage.io
import functools
import operator


from loguru import logger

logger.info('Import OK')

def magic(left, op, right):
   return op(left, right)

# define function to collect pixel location and intensity for a given mask, image combination
def pixel_collector(image_array, mask, mask_type=None, visualise=False, size=(1024, 1024), mask_cond=(operator.ne, 0)):
    """Obtains individual pixel coordinates and intensity values for an ROI.

    Parameters
    ----------
    image_array : 2D-array
        numpy array containing original image intensity values
    mask : 2D-array
        numpy array containing mask values for ROI of interest. Background pixels considered to be 0
    mask_type : [str], optional
        provided name of ROI type to be appended to coordinates df, by default None means column is not added
    visualise : bool, optional
        Determines whether to visualise the extracted ROI intensities via matplotlib, by default False
    size : tuple, optional
        If visualise=True, determines the x and y limits to recapitulate original image dimensions, by default (1024, 1024)
    mask_cond: tuple, optional
        Should take the form of (operation, value) e.g. to collect pixels where the mask is not equal to 0 (ne, 0) (default). Operations should be passed as per the standard operations module https://docs.python.org/3/library/operator.html.

    Returns
    -------
    DataFrame
        Pandas df containing x, y, intensity and optional mask_type columns
    """    

    pixel_array = np.where(magic(mask, *mask_cond), image_array, np.nan)
    # plt.imshow(pixel_array)
    coords = pd.DataFrame(pixel_array).unstack().reset_index().dropna()
    coords.columns = ['x', 'y', 'intensity']

    if mask_type != None:
        coords['mask_type'] = mask_type

    if visualise:
        # test visualisation, compare to plt.show
        fig, ax = plt.subplots(figsize=(20, 20))
        sns.scatterplot(coords['x'], coords['y'], hue=coords['intensity'], palette='magma_r', size=0.5,linewidth=0, alpha = 0.7)
        plt.ylim(size[0], 0)
        plt.xlim(0, size[1])

    return coords
