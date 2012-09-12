import numpy as np
from scipy.interpolate import UnivariateSpline

from .read_file import readCSVFile

def fit_length(csvdata, length='1e-4', k=1, s=0.):
    """ Fit a 1D spline from (x, y) csv extracted data.

    Retrieve the values it will be applied to from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'
    """
    length = float(length)
    if isinstance(csvdata, str):
        csvdata = readCSVFile(csvdata)
    x_name = csvdata.dtype.names[0]
    y_name = csvdata.dtype.names[1]

    print length, x_name, y_name, csvdata[x_name], csvdata[y_name]
    x = list(csvdata[x_name]/length)
    y = list(csvdata[y_name]/length)

    spline = UnivariateSpline(x, y, k=k, s=s)
    return spline



