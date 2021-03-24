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

    if not os.path.exists(f'{output_path}tifs/'):
        os.makedirs(f'{output_path}tifs/')

    folder_list = [sample for sample in os.listdir(input_path) if os.path.isdir(f'{input_path}{sample}')]

    for folder in folder_list:
        replicates = [sample for sample in os.listdir(f'{input_path}{folder}') if os.path.isdir(f'{input_path}{folder}/{sample}')]
        for replicate in replicates:
            replicate
            file_list = [filename for filename in os.listdir(f'{input_path}{folder}/{replicate}/') if '.tif' in filename]
            images = {}
            for filename in file_list:
                pathname = os.path.join(input_path, folder, replicate, filename)
                details = re.split(' ', filename)
                logger.info(details)
                if len(details) == 2:
                    new_name = f'{folder}_{replicate}_s1'+'.tif'
                    copyfile(pathname, output_path+'tifs/'+new_name)
                    images[1] = skimage.io.imread(f'{pathname}').transpose(1, 2, 0)
                    logger.info(f'{new_name}')
                elif len(details) == 3:
                    new_name = f'{folder}_{replicate}_s'+details[1].split('Pb')[1]+'.tif'
                    copyfile(pathname, output_path+'tifs/'+new_name)
                    images[int(details[1].split('Pb')[1])] = skimage.io.imread(f'{pathname}').transpose(1, 2, 0)
                    logger.info(f'{new_name}')
            # combine images into a single stack
            array_stack = np.dstack(tuple(images[key] for key in sorted(images.keys())))
            np.save(f'{output_path}{folder}_{replicate.strip("r")}.npy', array_stack)
                

jarvis(input_path='data/example_aggregate-FRAP/',
       output_path='results/aggregate-FRAP/initial_cleanup/')

