
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


def compute_K_from_laws(g, seminal_axial_conductivity_law, crown_axial_conductivity_law, lateral_axial_conductivity_law):
    """

    compute the axial conductance K versus the vertex position according to some laws and to the root types: crown, seminal,
    laterals

    :param g: (MTG)

    """
    K={}


    positions= g.property('relative_position')
    orders = g.property('order')

    for vid in g.vertices_iter(g.max_scale()):
            if g.label(vid)=='Crown':
                K[vid] = crown_axial_conductivity_law(positions[vid])
            elif g.label(vid) == 'Seminal' :
                K[vid] = seminal_axial_conductivity_law(positions[vid])
            elif g.label(vid) == 'Lateral' :
                    K[vid] = lateral_axial_conductivity_law(positions[vid])
            else: # for the collet and other unknown tpes
                K[vid] = seminal_axial_conductivity_law(positions[vid])

    g.properties()['K'] = K

    return g
