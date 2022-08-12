from math import pi
from collections import defaultdict

from openalea.mtg import *
#from openalea.mtg import algo

import numpy as np
from scipy.interpolate import UnivariateSpline
import pylab
from hydroroot.length import fit_law

def setting_k0_according_to_order(g, k0_pr, k0_lr):
    """
    set uniform radial conductivity to roots according to their order, to the primary root if oreder == 0,
    to the laterals otherwise
    :Parameters:
    	- g: (MTG) - the root architecture
    	- k0_pr: (float, microL/(s.MPa.m2)) - uniform radial conductivity of the primary root
    	- k0_lr: (float, microL/(s.MPa.m2)) - uniform radial conductivity of the laterals root

    :Returns:
        - g: (MTG) - the root architecture with k0_pr and k0_lr set
    """
    d = {k: k0_pr if g.property('order')[k] == 0 else k0_lr for k, v in list(g.property('order').items())}
    g.properties()['k0'] = d
    return g

def set_conductances(g, axial_pr, k0_pr, axial_lr = None, k0_lr = None):
    """
    Set the properties 'K_exp', 'K' and 'k', the axial conductance in [L^4 P^(-1) T^(-1)], the axial conductance in
    [L^3 P^(-1) T^(-1)] K_exp/segment_length and the radial conductivity
    if axial_lr is None, set 'K_exp' and 'K' whatever the roots order, otherwise set differently the root of order==1
    idem for the radial conductivity if k0_lr is not None

    :Parameters:
    	- g: (MTG)
    	- axial_pr: (list) - the axial conductance, list of 2 lists of floats
    	- k0_pr: (float) - the radial conductivity
    	- axial_lr:  (list) - if not None the axial conductance of the laterals, list of 2 lists of floats
    	- k0_lr: (float) - if not None the radial conductivity  of the laterals
    :Returns:
        - g (MTG)
    """
    xa, ya = axial_pr
    axial_conductivity_law = fit_law(xa, ya)

    # Compute K using axial conductance data
    if axial_lr is None:
        g = fit_property_from_spline(g, axial_conductivity_law, 'position', 'K_exp')
    else:
        K = {}
        xa, ya = axial_lr
        axial_conductivity_lr_law = fit_law(xa, ya)

        for v, k in list(g.property('order').items()):
            x = g.property('position')[v]
            if k == 0:
                K[v] = axial_conductivity_law(x)
            else:
                K[v] = axial_conductivity_lr_law(x)

        g.properties()['K_exp'] = K

    g = compute_K(g)  # Fabrice 2020-01-17: calculation of K in dimension [L^3 P^(-1) T^(-1)]

    if k0_lr is None: k0_lr = k0_pr

    g = setting_k0_according_to_order(g, k0_pr, k0_lr)

    g = compute_k(g, k0 = 'k0')

    return g

def compute_K_from_laws(g):
    K={}
    segment_length = 1e-4
    seminal_axial_conductivity_law = lambda x: 2135.05*x + 21338.63
    crown_axial_conductivity_law = lambda x: 9228.45*x + 42600.14
    lateral_axial_conductivity_law= lambda x: 2951.65*x + 2720.09

    positions= g.property('position')
    orders = g.property('order')

    for vid in g.vertices_iter(g.max_scale()):
            if g.label(vid)=='Crown':
                K[vid] = crown_axial_conductivity_law(positions[vid])
            elif g.label(vid) == 'Seminal' :
                if orders[vid] == 0:
                    K[vid] = seminal_axial_conductivity_law(positions[vid])
                else:
                    K[vid] = lateral_axial_conductivity_law(positions[vid])
            else:
                K[vid] = crown_axial_conductivity_law(positions[vid])

    g.properties()['K'] = K

    return g

def compute_K(g, scale_factor=1.):
    # Fabrice 2020-01-17: this calculation was done hydroroot.flux.run but that meant that the MTG was changed at each
    #                       flux calculation which is not relevant the MTG properties have to be fixed
    """
    Compute the conductance in dimension :math:`[L^3 P^{-1} T^{-1}]` from the experimental one which is in :math:`[L^4 P^{-1} T^{-1}]`

    :parameters:
        - g (MTG) : the root architecture
        - scale_factor (float) : a factor used for sensitivity analysis
    :return:
        - g (MTG) : with the property K calculated from property K_exp and the length of the vertex

    In each vertex K = K_exp * scale_factor / vertex_length
    """
    length = g.property('length')
    K_exp = g.property('K_exp')
    K = {}
    for vid in K_exp:
        K[vid] = K_exp[vid] / length[vid]
        K[vid] = K[vid] * scale_factor
    g.properties()['K'] = K
    return g



def poiseuille(radius, length, viscosity=1e-3):  # DEPRECATED
    """
    Deprecated
    Compute a conductance of a xylem element based on their radius and length.
    
    :parameters:
        - radius (float) : (m) radius of a xylem tube

    length (float) : (m) length of a xylem element

    viscosity (float) : (Pa.s) dynamic viscosity of the liquid
    
    The poiseuille formula is:
        :math:` conductance = \frac{\pi r^4}{8 \mu L }` 
        with :math:`r` the radius of a pipe, 
        :math:`\mu` the viscosity of the liquid,
        :math:`L` the length of the pipe.
        
    .. seealso:: http://en.wikipedia.org/wiki/Poiseuille
    """
    return pi*(radius**4) / ( 8 * viscosity * length)


def compute_k(g, k0 = 300.):
    """
    Compute the radial conductances (k) of each segment of the MTG.

    :parameters:
        - g (MTG) : the RSA
        - k0 (float) or (string): the radial conductance for one element of surface in :math:`micronL/s.MPa.m^2`

    if k0 == 'k0': calculation using k0 from g.property('k0')
    if k0 is a float: value used in calculation

    k radial (microL/s.MPa) calculation:
        k = 2 \pi r l k0 with l and r the segment length and radius

    """
    #print 'entering radial k fitting'

    radius = g.property('radius')
    length = g.property('length')
    kr={}
    if k0 == 'k0':
        k0 = g.property('k0')
        kr = dict((vid, radius[vid] * 2 * pi * length[vid] * k0[vid]) for vid in g.vertices(scale=g.max_scale()))
    else:
        kr = dict((vid, radius[vid] * 2 * pi * length[vid] * k0) for vid in g.vertices(scale=g.max_scale()))

    g.properties()['k'] = kr
    #print 'exiting radial k fitting'
    return g


def compute_K_from_Poiseuille(g, nb_xylem=5, radius_scale = 1/10.):  # DEPRECATED
    # Fabrice 2020-01-17: changed the function name from "compute_K" to "compute_K_from_Poiseuille" because this name
    #                     is now used to calculate the real conductance in [L^3 P^(-1) T^(-1)] from the experimental one
    #                     in [L^4 P^(-1) T^(-1)]
    """
    Deprecated
    Compute the axial conductances (K) in a MTG according to Poiseuille law.

    The conductance depends on the radius of each xylem pipe, the number of xylem pipes,
    and on the length of a root segment.

    radius_scale allows to compute the radius of a xylem pipe from the radius of a root segment.
    """

    radius = g.property('radius_xylem')
    if not radius:
        full_radius = g.property('radius')
        radius = dict( (vid,r*radius_scale) for vid,r in full_radius.items())
    nb_xylem = g.property('nb_xylem')
    length= g.property('length')
    if not nb_xylem:
        nb_xylem = defaultdict(lambda : 5)
    K = dict((vid, nb_xylem[vid]*poiseuille(radius[vid], length[vid])) 
                for vid in g.vertices(scale=g.max_scale()))
    g.properties()['K'] = K
    return g

def fit_property(g, x, y, prop_in, prop_out, s=3.): 
    """
    Deprecated
    Fit a 1D spline from x, y data.

    Retrieve the values from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'
    """

    spline = UnivariateSpline(x, y, s=s)
    keys = list(g.property(prop_in).keys())
    x_values = np.array(list(g.property(prop_in).values()))

    y_values = spline(x_values)

    g.properties()[prop_out] = dict(list(zip(keys,y_values)))

    xx = np.linspace(0,1,1000)
    yy = spline(xx)

    pylab.clf()
    pylab.plot(x, y)
    pylab.plot(xx, yy)
    pylab.show()

    #print 'Update figure ', yy.min(), yy.max()
    return g


def fit_property_from_spline(g, spline, prop_in, prop_out): 
    """
    compute a property from another one using a spline transformation.

    Retrieve the values from the prop_in of the MTG.
    And evaluate the spline to compute the property prop_out

    :parameters:
        - g (MTG)
        - spline (class scipy.interpolate.UnivariateSpline) : 1-D smoothing spline fit to a given set of data points
        - prop_in (string) : the property data
    """

    #spline = UnivariateSpline(x, y, s=s)
    keys = list(g.property(prop_in).keys())
    x_values = np.array(list(g.property(prop_in).values()))

    y_values = spline(x_values)

    g.properties()[prop_out] = dict(list(zip(keys, y_values)))

    return g


def fit_property_from_csv(g, csvdata, prop_in, prop_out, k=1., s=0., plot=False, direct_input=None):
    """
    Deprecated
    Fit a 1D spline from (x, y) csv extracted data or from direct input dictionnary

    Retrieve the values it will be applied to from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'

    Toggle plot option to visualize the spline fit
    """
    #print 'entering K fitting'    

    from .read_file import readCSVFile
    if isinstance(csvdata, str):
        csvdata = readCSVFile(csvdata)

    if direct_input is None:
        x_name = csvdata.dtype.names[0]
        y_name = csvdata.dtype.names[1]
        x = list(csvdata[x_name])
        y = list(csvdata[y_name])
    else:
        x, y = [], []
        for key in sorted(direct_input.keys()):   # dictionnary key are not ordered by default
            x.append(key)
            y.append(direct_input[key])

    spline = UnivariateSpline(x, y, k=k, s=s)
    fit_property_from_spline(g, spline, prop_in, prop_out)

    if plot:
        x_n = np.array(list(g.property(prop_in).values()))
        #y_n = np.array(g.property(prop_out).values())

        xx = np.linspace(min(x_n), max(x_n), 1000)
        yy = spline(xx)

        # plot the reference (x_values,y_values) data and the fitted spline
        pylab.clf()
        pylab.plot(x, y, 'x')
        #pylab.plot(x_values, y_values)
        pylab.plot(xx, yy, '-')
        pylab.show()

        #print 'Update figure ', xx.min(), yy.max()

    #print 'exiting K fitting'

    return g


def fit_K(g, s=0.):   # DEPRECATED
    """
    Deprecated
    """
    x = np.linspace(0.,1.,100)
    y = np.linspace(50, 500, 100)+100*np.random.random(100)-50

    if s == 0.:
        s = None
    fit_property(g,x,y,'relative_position', 'K', s=s)


    return g


# Below these function does not do very complicated things but are used in most of the scripts

def radial(v = 92, acol = [], scale = 1):
    """
    create a list of uniform value v*scale of the same length than acol given in arguments
    the purpose is to return  x-y data in a form of two lists
    called radial because used to get x-y data with uniform y (radial conductivity) values

    :Parameters:
    	- v: (float)
    	- acol: (float list)
    	- scale: (float) -
    :Returns:
        - xr, yr (list)
    """
    xr = acol[0]  # at this stage kr constant so the same x than Ka
    yr = [v * scale] * len(xr)

    return xr, yr


def radial_step(kmin = 92, factor = 1.0, x_step = None, x_base = 1., dx = 1.0e-4, scale = 1.0):
    """
    Radial k as a step function.
    Maximum value from the tip x=0 to the x_step, then minimum value from x_step + dx to x_base.
    The maximum value = kmin * factor

    :Parameters:
        - kmin: (Float), the minimum value in microL/(s.MPa.m**2)
        - factor: (Float), see above
        - x_step: (Float), the distance from tip where the step is in meter
        - x_base: (Float), the maximum distance from tip in meter
        - dx: (Float), elementary distance in m
        - scale: (Float), a scale facteur

    :Returns:
        - xr: (list), list of distance from the tip
        - xr: (list), list of radial k values
    """

    xr = [0.0]
    yr = [kmin * factor * scale]
    if x_step is not None:
        xr.append(x_step)
        yr.append(kmin * factor * scale)
        xr.append(x_step + dx)
        yr.append(kmin * scale)
    xr.append(x_base)
    yr.append(kmin * scale)

    return xr, yr


def axial(acol = [], scale = 1):
    """
    the purpose is to give in arguments a set of 2 lists representing a x-y data and to return it with y*scale

    :Parameters:
    	- acol: (list) - list of two float list of the same length
    	- scale: (float) - the number that multiply the 2d list
    :Returns:
        - x, y (list) - 2 lists
    """
    x, y = acol
    y = [a * scale for a in y]
    return x, y
