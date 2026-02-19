
""""
# Conductance module for HydroRoot.Millet

"""

import numpy as np
from scipy.interpolate import UnivariateSpline

def fit_property_from_spline(g, spline, prop_in, prop_out): 
    """ compute a property from another one using a spline transformation.

    Retrieve the values from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'
    """

    #spline = UnivariateSpline(x, y, s=s)
    keys = g.property(prop_in).keys()
    x_values = np.array(list(g.property(prop_in).values())) # g.property(prop_in).values() return a dict_values object, need to be converted a list first

    y_values = spline(x_values)

    g.properties()[prop_out] = dict(zip(keys, y_values))

    return g

def compute_K_from_laws(g, seminal_law=None,crown_law=None, lateral_law=None,):
    """

    Compute axial hydraulic conductance (K) for each MTG vertex from
    user-defined laws as a function of the vertex position along the axis.


    Parameters
    ----------
    g : MTG
        Multi-scale Tree Graph containing root architecture.

    seminal_law : callable, optional
        Function f(x) -> K defining axial conductance for seminal roots     

    crown_law : callable, optional
        Function f(x) -> K defining axial conductance for crown roots.

    lateral_law : callable, optional
        Function f(x) -> K defining axial conductance for lateral roots

    Returns
    -------
    g : MTG
        The input MTG with a new vertex property 'K' added.

    """
    K={}

    positions= g.property('relative_position')
    orders = g.property('order')

    for vid in g.vertices_iter(g.max_scale()):
        label = g.label(vid)
        pos = positions[vid]
        if label == "Crown":
            K_val = crown_law(pos)

        elif label == "Seminal":
            K_val = seminal_law(pos)

        elif label == "Lateral":
            K_val = lateral_law(pos)

        else:
            # fallback (e.g., collet or unknown types)
            K_val = seminal_law(pos)

        K[vid] = K_val

    # Store property in the MTG
    g.properties()['K'] = K

    return g
