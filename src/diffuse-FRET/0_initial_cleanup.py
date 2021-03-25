import os
import re
from shutil import copyfile
import numpy as np
import skimage.io


from loguru import logger
logger.info('Import ok')

def jarvis(input_path, output_path):
    """Simple filename cleansing for exported TIF filenames (as these are normally a mess following export via ImageJ)

    Parameters
    ----------
    input_path : str
        folder where input files are stored
    output_path : str
        folder where renamed files should be copied to.
    """    
    
    if not os.path.exists(f'{output_path}'):
        os.makedirs(f'{output_path}')

    folder_list = [sample for sample in os.listdir(input_path) if os.path.isdir(f'{input_path}{sample}')]

    for folder in folder_list:

        file_list = [filename for filename in os.listdir(f'{input_path}{folder}/') if '.tif' in filename]
        mutant = '_'.join(folder.split(' '))

        for x, filename in enumerate(file_list):
            pathname = os.path.join(input_path, folder, filename)
            new_name = f'{output_path}{mutant}_{x}.tif'
            copyfile(pathname, new_name)
            # array_stack = skimage.io.imread(f'{pathname}').transpose(1, 2, 0)
            logger.info(f'{new_name}')
                

if __name__ == "__main__":

    jarvis(input_path='data/example_diffuse-FRET/',
           output_path='results/example_diffuse-FRET/initial_cleanup/')
