# -*- coding: utf-8 -*-
"""
Define a length law of LRs using equal amplitudes method from measured data.

Created on Tue Feb 17 12:56:51 2015

@author: ndour
"""

import random
import numpy as np
import pandas as pd
import pylab

from scipy.optimize import curve_fit

from openalea.mtg import algo
from openalea.mtg.algo import orders, axis
from openalea.mtg.traversal import pre_order2_with_filter, post_order2

from openalea.hydroroot import length, radius
from openalea.hydroroot import conductance



def read_data_deprecated(data_xy, scale_x, scale_y,sort = True):
    # sort by position
    #  8/27/25 xy.sort(axis=0) does sort both x and y we do not want that
    xy = data_xy
    if sort == True:
        xy.sort(axis=0)
    xy = xy.tolist()
    x, y = zip(*xy)  # Separate x and y coordiantes
    x = np.array(x) * scale_x
    # x = np.array(x)
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    return x,y

def read_data(data_xy, scale_x, scale_y,sort = True):
    # sort by position

    if sort == True:
        data_xy = data_xy.sort_values(data_xy.columns[0])

    x = np.array(data_xy[data_xy.columns[0]]) * scale_x
    y = np.array(data_xy[data_xy.columns[1]]) * scale_y
    x = x.tolist()
    y = y.tolist()

    return x,y


def expovariate_law(data_xy, size=5e-2, segment_length =1e-4, plot=False):
    """
    Fit a spline law from measured data by adding stochasticity.

    To compute the law, data are first sampled using equal amplitude method. Then, a mean is computed on each sample.
    Then, an expovariate distribution is simulated from the mean of each sample.
    Finally a spline is interpolated based on the simulated data.

    Parameters
    ==========
        - data_xy: DatFrame
        - size:
        - scale_x:
        - scale_y:
        - plot:

    Example
    =======

        >>> filename = 'lr_length_law_data.csv'
        >>> xy = read_CSVFile(filename)
        >>> law = fit_law(data_xy=xy, size=5e-2)
    """
    #return x and y as lists
    x,y = read_data(data_xy, scale_x=1e-2, scale_y=1e-2)


    X, values = discretize(x, y, size=size)

    Y = [np.mean(ys) for ys in values]
    YY = [(random.expovariate(1. / v) if v > 0 else 0.) for v in Y]
  

    if plot:
        Y_max = [max(ys) for ys in values]
        Y_min = [min(ys) for ys in values]

        pylab.plot(x, y)
        pylab.plot(X, Y_max)
        pylab.plot(X, Y_min)
        pylab.plot(X, YY)

    #Convert length(m) to number of vertex
    X = list(np.array(X) / segment_length)
    YY = list(np.array(YY) / segment_length) 
    

    law = length.fit_law(X, YY)

    return law


def discretize(x, y, size=5e-2):
    """ Discretize by 5cm-intervals by using equal amplitudes method"""

    m, M = min(x), max(x)

    nb_class = int((M - m) / size)

    points = [(m + i * size) for i in range(1, nb_class - 1)]

    points.insert(0, m)
    points.append(M)
    nb_points = len(points)
    intervals = [(points[i], points[i + 1]) for i in range(nb_points) if i != (nb_points - 1)]

    ys = [[y[i] for i, p in enumerate(x) if p1 <= p <= p2] for p1, p2 in intervals]

    ys.insert(0, [0.])

    return points, ys

######################################################################################################



def diameter_law(data_xy, function=None, segment_length = 1e-4, plot= True):

    """this function checks the the best paramters of the fonction that fit the data"""

    x,y = read_data(data_xy, scale_x=1e-2, scale_y=1e-6)
    #Convert x to vertex
    
    x = list(np.array(x)/segment_length)
    #The popt argument are the best-fit paramters for a and b
    popt, pcov = curve_fit(function, x, y)
    #Simulate y-data using best-fit parameters
    yy= function(x,*popt)


    if plot:
        pylab.plot(x,y,'ko', label="Original Data")
        pylab.plot(x, yy, 'r-', label="Fitted Curve")
        pylab.legend()

    def f(k,fn = function):
        return fn(k,*popt)

    law = f
    return law

def diameter_law2(data_xy, function=None, plot= True):

    """this function checks the the best paramters of the fonction that fit the data"""

    x,y = read_data(data_xy, scale_x=1e-2, scale_y=1e-6)
    #convert to float
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Normalize x to have a relative position
    #x = list(np.array(x)/np.max(np.array(x)))
    x = x / np.max(x)
    #The popt argument are the best-fit paramters for a and b
    popt, pcov = curve_fit(function, x, y)
    #Simulate y-data using best-fit parameters
    yy= function(x,*popt)


    if plot:
        pylab.plot(x,y,'ko', label="Original Data")
        pylab.plot(x, yy, 'r-', label="Fitted Curve")
        pylab.legend()

    def f(k,fn = function):
        return fn(k,*popt)

    law = f
    return law



def compute_radius_from_laws(g,
                     seminal_RootDiameter_law=None,
                     crown_RootDiameter_law=None,
                     lr_RootDiameter_law=None,
                     segment_length = 1e-4):
    
    
    #g= radius.compute_length(g)
    #g= radius.compute_relative_position(g)
    #position = g.property('position')
    #order= g.property('order')

    for v in g.vertices_iter(g.max_scale()):
        node = g.node(v)
        pos = node.position/segment_length #convert position which is expressed in m to number of vertex
        #We use a logarithmic law that takes position as input to estimate diameter so the position can be null
        if pos <= 1 or node._vid in [0,1,2,3]: #pos <= 1 to avoid ln(0) in RootDiameter_law
            node.radius = 270.*1e-6
        else:
            if node.label == 'Crown':
                node.radius = crown_RootDiameter_law(pos)
            elif node.label == 'Seminal':
                if node.order == 0:
                    node.radius = seminal_RootDiameter_law(pos)
                else:
                    node.radius = lr_RootDiameter_law(pos)

    return g



# ************** Here we define evolution diameter laws of different root types************************


def compute_diameter_from_laws(g,
                                seminal_RootDiameter_law=None,
                                crown_RootDiameter_law=None,
                                lr_RootDiameter_law=None
                               ):
    """
    calculate and add diameter as a propoerty to the MTG

    :params:
    - g: MTG
    - seminal_RootDiameter_law: law of seminal root diameter
    - crown_RootDiameter_law: law of crown root diameter
    - lateral_RootDiameter_law: law of lateral root diameter

    :returns:
    g: MTg with diameter property
    """

    diameters = {}
    order= g.property('order')
    positions = g.property('relative_position')
    for vid in g.vertices_iter(g.max_scale()):
        # collet or Seminal
        if order[vid] == 0:
            diameters[vid] = seminal_RootDiameter_law(positions[vid])
        else:
            if g.label(vid)== "Crown":
                diameters[vid] = crown_RootDiameter_law(positions[vid])
            else: # laterals
                diameters[vid] = lr_RootDiameter_law(positions[vid])  

    g.properties()['diam'] = diameters

    return g


def compute_radius_from_diameter(g):
    """
    Calculate and add radius to MTG as a property from diameter values
    """
    radius = {}
    diameters= g.property('diam')
    for i,d in diameters.items():
        radius[i] = d/2.
    g.properties()['radius'] = radius

    return g

def compute_diameter(g,
                     seminal_RootDiameter_law=None,
                     crown_RootDiameter_law=None,
                     lr_RootDiameter_law=None,
                     segment_length = 1e-4):


    diameters = {}
    order= g.property('order')
    g = radius.compute_length(g, segment_length)
    g= radius.compute_relative_position(g)
    positions = g.property('relative_position')
    for vid in g.vertices_iter(g.max_scale()):
            if g.label(vid)=='Crown':
                if positions[vid] < 1:
                    diameters[vid] = 380.0
                else:
                    diameters[vid] = crown_RootDiameter_law(positions[vid])
            elif g.label(vid) == 'Seminal':
                if order[vid] == 0:
                    if positions[vid] < 1:
                        diameters[vid] = 260.0
                    else:
                        diameters[vid]= seminal_RootDiameter_law(positions[vid])
                else:
                    if positions[vid] < 1:
                        diameters[vid]= 240.0
                    else:
                        diameters[vid] = lr_RootDiameter_law(positions[vid])
            else : #for the collet
                diameters[vid] = 380.0


    g.properties()['diam'] = diameters

    return g

###################################################################

def radius_from_computed_diameters(g):
    radius = {}
    diameters= g.property('diam')
    for i,j in diameters.items():
        radius[i] = (j/2.)*1e-4
    g.properties()['radius'] = radius

    return g



###############################################################################

def age_law(data_xy, segment_length = 1e-4,scale_x =1, scale_y = 1.e-2, plot= True):
    """compute age property from position of each vertex """

    x,y = read_data(data_xy, scale_x=scale_x, scale_y=scale_y, sort=False)
    #Convert x to vertex
    
    y = list(np.array(y)/segment_length)
    x= list(x)
    #The popt argument are the best-fit paramters for a and b
    #popt, pcov = curve_fit(function, x, y)

    law = length.fit_law(y,x)

    xx = law(y)


    if plot:
        pylab.plot(x,y,'ko', label="Original Data")
        pylab.plot(xx, y, 'r-', label="Fitted Curve")
        pylab.legend()

    return law


################################################################################

def developmental_age(g, nude_tip_length=15):
    """ Compute the developmental age of each vertex.

    A new property is added to the MTG that represent the age of apparition of the segment.
    It allows to express a dynamic as a parametrisation.
    """
    age = {}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    is_seminal = lambda v: g.label(v) != 'Crown'

    for v in post_order2(g, 1, pre_order_filter=is_seminal):
        if g.is_leaf(v):
            age[v] = 0  
        else:
            age[v] = max(age.get(cid,0) for cid in g.children(v))+1

    age_max = age[1]

    k, v = np.array(list(age.keys())), np.array(list(age.values()))
    v = age_max - v

    age = dict(zip(list(k), list(v)))

    delta = nude_tip_length

    for v in pre_order2_with_filter(g, 1, pre_order_filter=is_seminal):
        pid = g.parent(v)
        if g.edge_type(v) =='+':
            age_p, age_v = age[pid], age[v]
            if age_p +delta < age_v:
                age[v] = age_p + delta  # replace by a stat distribution rather than min
        elif pid is not None:
            age[v] = age[pid] + 1

    # Crown roots : Same for seminal with a delta
    # and a different speed?
    delta_crown_age = 2344

    root_crowns = [v for v in g if g.label(v) == 'Crown' and g.label(v) != g.label(g.parent(v))]

    speed = 1.
    for v in root_crowns:
        date = delta_crown_age
        for cid in axis(g, v):
            age[cid] = date
            date += speed

    g.properties()['age'] = age
    return g



def developmental_age2(g, nude_tip_length=0.0206):
    """ Compute the developmental age of each vertex.

    A new property is added to the MTG that represent the age of apparition of the segment.
    It allows to express a dynamic as a parametrisation.
    """
    age = {}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            # convert nude_tip_length to nb of vertices: for a 75cm root, the length of nude tip is 15cm
    nude_tip_length = (nude_tip_length * 0.15) / 0.75

    #--------------------------------------
    # Seminal
    #--------------------------------------
    is_primary = lambda v: g.label(v) == 'Seminal' or g.label(v)=='collet'

    for v in post_order2(g, 1, pre_order_filter=is_primary):
        if g.is_leaf(v):
            age[v] = 0  
        else:
            age[v] = max(age.get(cid,0) for cid in g.children(v))+1

    age_max = age[1]

    k, v = np.array(list(age.keys())), np.array(list(age.values()))
    v = age_max - v

    age = dict(zip(list(k), list(v)))

    delta = nude_tip_length

    for v in pre_order2_with_filter(g, 1, pre_order_filter=is_primary):
        pid = g.parent(v)
        if g.edge_type(v) =='+':
            age_p, age_v = age[pid], age[v]
            if age_p +delta < age_v:
                age[v] = age_p + delta  # replace by a stat distribution rather than min
        elif pid is not None:
            age[v] = age[pid] + 1

    #---------------------------------------------------------------
    # Crown: they appear 6 days after germination from the root base
    #---------------------------------------------------------------

    delta_crown_age = 6

    root_crowns = [v for v in g if g.label(v) == 'Crown' and g.label(v) != g.label(g.parent(v))]

    speed = 1.
    for v in root_crowns:
        date = delta_crown_age
        for cid in axis(g, v):
            age[cid] = date
            date += speed


    #--------------------------------------------------------------------
    # Lateral: they appear 6 days after germination from the primary root
    #--------------------------------------------------------------------

    delta_crown_age = 6

    root_laterals = [v for v in g if g.label(v) == 'Lateral' and g.label(v) != g.label(g.parent(v))]

    speed = 1.
    for v in root_laterals:
        date = delta_crown_age
        for cid in axis(g, v):
            age[cid] = date
            date += speed



    # Store  the property in the MTG
    g.properties()['age'] = age
    return g
###########################################################################################################


def add_soil_old(g, soil_data, segment_length = 1.e-4):

    """ add a soil a hetergeneous water potential """

    x,y = soil_data

    # Compute absolute z coordinate and normalize
    vids = g.property('xyz').keys()
    zs = np.array([pt.z * segment_length for pt in g.property('xyz').values()])
    zs/=zs.max()
    zs = zs.tolist()
    # Fit data on z coordinate to compute psi_e on each vertex
    g.properties()['height'] = dict(zip(vids,zs))
    soil_law = length.fit_law(x,y)
    g = conductance.fit_property_from_spline(g, soil_law, 'height', 'psi_e')

    return g


def add_soil(g, soil_data, segment_length = 1.e-4):

    """ add a soil a hetergeneous water potential """

    x,y = soil_data

    # Compute absolute z coordinate and normalize
    vids = g.property('xyz').keys()
    zs = np.array([np.abs(pt.z) * segment_length for pt in g.property('xyz').values()])
    zs/=zs.max()
    zs = zs.tolist()
    # Fit data on z coordinate to compute psi_e on each vertex
    g.properties()['height'] = dict(zip(vids,zs))
    soil_law = length.fit_law(x,y)
    g = conductance.fit_property_from_spline(g, soil_law, 'height', 'psi_e')

    return g