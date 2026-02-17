
""""
# Conductance module for HydroRoot.Millet

"""

import numpy as np
#from scipy.interpolate import UnivariateSpline

def fit_property_from_spline(g, spline, prop_in, prop_out): 
    """ compute a property from another one using a spline transformation.

    Retrieve the values from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'
    """

    #spline = UnivariateSpline(x, y, s=s)
    keys = g.property(prop_in).keys()
    x_values = np.array(g.property(prop_in).values())

    y_values = spline(x_values)

    g.properties()[prop_out] = dict(zip(keys, y_values))

    return g

