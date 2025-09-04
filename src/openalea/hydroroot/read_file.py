import numpy as np
import pandas as pd

def readCSVFile(filename):
    """Read and extract data from a csv file, supposed that the data is stored in 2 columns.

    :param filename: string
    :returns: - data (array) - record array of (x, y) values, column headers recorded in dtype

    """
    data = np.recfromcsv(filename,delimiter=';')
    return data

def read_archi_data(fn):
    """Read a csv (tab separated) file with the architecture in the following format

    :param fn: (string) - the architecture filename in csv format with tab as delimiter
    :returns:
        - DataFrame, with columns from the csv plus:
            - a column 'db', the distance to base in m, and
            - a column 'lr', the lateral root lengths in m
            - a column 'radius', if 'averaged_diameter_(mm)' the radius in m

    :file structure:

    +--------------------------+----------------------------+-------+-------------------------+
    |'distance_from_base_(mm)' | 'lateral_root_length_(mm)' | order | 'averaged_diameter_(mm)'|
    +--------------------------+----------------------------+-------+-------------------------+
    |float                     | float                      | string| float                   |
    +--------------------------+----------------------------+-------+-------------------------+

    - order = 1 for laterals of 1st order ob the primary
    - order = n-m for the lateral number m on the lateral number n of 1st order
    - order = n-m-o for the lateral number o of the previous one

    Column 'averaged_diameter_(mm)' is optional

    Each branch finish with a nude part, i.e. a distance from base (the tip value) and a zero length
    The 'averaged_diameter_(mm)' is not mandatory, if present allow calculation of radii for each laterals
    if not the radii of PR and LR will be calculated in routine radius.ordered_radius.
    the resulting dataframe must be in meter and have the following column names: 'db', 'lr', 'order' and 'radius' (optional)
    these names are used in mtg_from_aqua_data to build the MTG

    """
    df = pd.read_csv(fn, sep = '\t', dtype = {'order': str})
    df['db'] = df['distance_from_base_(mm)'] * 1.e-3
    df['lr'] = df['lateral_root_length_(mm)'] * 1.e-3

    if 'averaged_diameter_(mm)' in df:
        df['radius'] = df['averaged_diameter_(mm)'] * 0.5e-3

    return df
