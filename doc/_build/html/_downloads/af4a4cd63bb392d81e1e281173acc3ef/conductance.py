from math import pi
from collections import defaultdict

from openalea.mtg import *
#from openalea.mtg import algo

import numpy as np
from scipy.interpolate import UnivariateSpline
import pylab



def compute_K_from_laws(g):
    """
    Deprecated
    Compute the axial conductance from empirical laws (hard coded) according to the distance from the tip and to the
    root type, either seminal, crown or lateral

    :parameters:
        - g (MTG)

    set property K
    """
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
    #                     is now used to calculate the real conductance in [L^3 P^{-1} T^{-1}] from the experimental one
    #                     in [L^4 P^{-1} T^{-1}]
    """ Compute the axial conductances (K) in a MTG according to Poiseuille law.

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
    """ Fit a 1D spline from x, y data.

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
    """ compute a property from another one using a spline transformation.

    Retrieve the values from the prop_in of the MTG.
    And evaluate the spline to compute the property 'prop_out'
    """

    #spline = UnivariateSpline(x, y, s=s)
    keys = list(g.property(prop_in).keys())
    x_values = np.array(list(g.property(prop_in).values()))

    y_values = spline(x_values)

    g.properties()[prop_out] = dict(list(zip(keys, y_values)))

    return g


def fit_property_from_csv(g, csvdata, prop_in, prop_out, k=1., s=0., plot=False, direct_input=None):
    """ Fit a 1D spline from (x, y) csv extracted data or from direct input dictionnary

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
    x = np.linspace(0.,1.,100)
    y = np.linspace(50, 500, 100)+100*np.random.random(100)-50

    if s == 0.:
        s = None
    fit_property(g,x,y,'relative_position', 'K', s=s)


    return g
