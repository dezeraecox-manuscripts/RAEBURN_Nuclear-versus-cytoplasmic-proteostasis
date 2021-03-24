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

    file_list = [filename for filename in os.listdir(f'{input_path}') if '.tif' in filename]

    for filename in file_list:
        pathname = os.path.join(input_path, filename)
        new_name = f"{output_path}{filename.replace('.lif - ', '_').replace('_5x-', '_')}"
        copyfile(pathname, new_name)
        logger.info(f'{new_name}')
                
if __name__ == "__main__":

    jarvis(input_path='data/example_chaperone-localisation/',
        output_path='results/chaperone_localisation/initial_cleanup/')

