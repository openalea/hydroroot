import numpy as np
from scipy.interpolate import UnivariateSpline

from .read_file import readCSVFile


def fit_length(csvdata, length='1e-4', k=1, s=0.):
    """Fit a 1D spline from (x, y) csv extracted data.
    
    Retrieve the values it will be applied to from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'

    :param csvdata: 
    :param length: number dividing data if > 0 (Default value = '1e-4')
    :param k: (int) - degree of the smoothing spline (Default value = 1)
    :param s: (float) - positive smoothing factor used to choose the number of knots (Default value = 0)
    :return:
        - spline object

    """
    length = float(length)
    if isinstance(csvdata, str):
        csvdata = readCSVFile(csvdata)
    x_name = csvdata.dtype.names[0]
    y_name = csvdata.dtype.names[1]

    return fit_law(csvdata[x_name], csvdata[y_name], scale=length, k=k, s=s)


def fit_law(x, y, scale=0., k=1, s=0, **kwds):
    """
    Return a spline interpolation of y(x)

    :param x: (float list)
    :param y:  (float list)
    :param scale: (float) - number dividing x and y if > 0 (Default value = 0.)
    :param k: (int) - degree of the smoothing spline (Default value = 1)
    :param s: (float) - positive smoothing factor used to choose the number of knots (Default value = 0)
    :param kwds: additional arguments see scipy doc of UnivariateSpline
    :return:
        - spline object
    """
    if scale:
        x = list(np.array(x) / scale)
        y = list(np.array(y) / scale)

        #print "DEBUG: ", scale, x, y
    spline = UnivariateSpline(x, y, k=k, s=s, **kwds)
    return spline

def diff(law1, ref_law):
    """
    deprecated
    Calculate the difference between the inetgrale of the two laws
    :param law1: scipy spline object
    :param ref_law: scipy spline object

    """
    knots = law1.get_knots()

    interval_def = (knots[0], knots[-1])
    integral1 = law1.integral(*interval_def)
    integral_ref = ref_law.integral(*interval_def)

    return integral1-integral_ref
