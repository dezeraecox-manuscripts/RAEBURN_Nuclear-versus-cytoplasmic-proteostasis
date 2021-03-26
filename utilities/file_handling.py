import os, re
import pandas as pd
import numpy as np

from loguru import logger

logger.info('Import OK')


def df_to_excel(output_path, sheetnames, data_frames):
    """Saves list of dataframes to a single excel (xlsx) file.
    Parameters
    ----------
    output_path : str
        Full path to which xlsx file will be saved.
    sheetnames : list of str
        descriptive list of dataframe content, used to label sheets in xlsx file.
    data_frames : list of DataFrames
        DataFrames to be saved. List order must match order of names provided in sheetname.
    Returns
    -------
    None.
    """
    if not output_path.endswith('.xlsx'):
        output_path = output_path+'Results.xlsx'
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    # Convert the dataframe to an XlsxWriter Excel object.
    for x in range(0, len(sheetnames)):
        sheetname = sheetnames[x]
        data_frame = data_frames[x]
        data_frame.to_excel(writer, sheet_name=sheetname)
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

