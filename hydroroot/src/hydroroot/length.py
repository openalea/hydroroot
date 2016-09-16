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

    return fit_law(csvdata[x_name], csvdata[y_name], scale=length, k=k, s=s)


def fit_law(x, y, scale=0., k=1, s=0.):
    if scale:
        x = list(np.array(x) / scale)
        y = list(np.array(y) / scale)

        print "DEBUG: ", scale, x, y
    spline = UnivariateSpline(x, y, k=k, s=s)
    return spline
