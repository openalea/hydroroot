###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to compute architecture or to construct it from texte file
#   Display them with a colormap based on a MTG property ('j', 'order', etc.)
#       argument --prop passed through command line
#   If order chosen then display root with a color according to there order
###############################################################################

######
# Imports

# VERSION = 2
import copy
from binascii import hexlify as _hexlify # F. Bauget 2021-12-14  # "from random import _hexlify, _urandom"
from os import urandom as _urandom # F. Bauget 2021-12-14  # "from random import _hexlify, _urandom"

import pandas as pd

from openalea.mtg import traversal

from hydroroot import radius, markov
from hydroroot.law import histo_relative_law
from hydroroot.generator.measured_root import mtg_from_aqua_data

# read architecture file
def read_archi_data(fn):
    """
    Read a csv (tab separated) file with the architecture in the following format
        |'distance_from_base_(mm)' | 'lateral_root_length_(mm)' | order |
        |float | float | string|
        order = 1 for laterals of 1st order ob the primary
        order = n-m for the lateral number m on the lateral number n of 1st order
        order = n-m-o for the lateral number o of the previous one
        etc.
        Each branch finish with a nude part, i.e. a distance from base (the tip value) and a zero length

    :param fn: string - the architecture filename in csv format

    :return: DataFrame
    """
    df = pd.read_csv(fn, sep = '\t', dtype = {'order': str})
    df['db'] = df['distance_from_base_(mm)'] * 1.e-3
    df['lr'] = df['lateral_root_length_(mm)'] * 1.e-3

    return df

#################################################################################
# MTG construction either from data reconstructed() or generated from parameters
#   generated()
#################################################################################

def generate_g(seed = None, length_data = None, branching_variability = 0.25,
               delta = 2e-3, nude_length = 2e-3, primary_length = 0.13, segment_length = 1e-4, order_max = 4):
    """generate a MTG according to the input parameters
    Author: F. Bauget, based on routine writtent by C. Pradal, Y. Boursiac
    Date: 2019-12-10

    Parameters
    ==========
        - seed: the seed for the random generator in the markof chain
        - length_data: pandas dataframe columns names 'LR_length_mm', 'relative_distance_to_tip' sorted by 'relative_distance_to_tip'
        - branching_variability: probability of ramification at exact mean branching position
        - branching_delay: reference distance between successive branching axis
        - nude_length: length at root tip with no ramification
        - primary_length: primary root length
        - segment_length: length of the vertices, default 1.e-4
        - order_max: maximum lateral roots order

    Returns
    =======
        - g: MTG with the following properties set: edge_type, label, position
    """

    # nude length and branching delay in terms of number of vertices
    nb_nude_vertices = int(nude_length / segment_length)
    branching_delay = int(delta / segment_length)

    nb_vertices = int(primary_length / segment_length)

    length_max_secondary = length_data[0].LR_length_mm.max() * 1e-3  # in m

    law_order1 = length_law(length_data[0], scale_x = primary_length / 100., scale = segment_length)
    law_order2 = length_law(length_data[1], scale_x = length_max_secondary / 100., scale = segment_length)

    g = markov.markov_binary_tree(
        nb_vertices = nb_vertices,
        branching_variability = branching_variability,
        branching_delay = branching_delay,
        length_law = [law_order1, law_order2],
        nude_tip_length = nb_nude_vertices,
        order_max = order_max,
        seed = seed)
    return g

###############################################################################
# data
###############################################################################
def length_law(pd, scale_x = 1 / 100., scale_y = 1., scale = 1e-4, uniform = 'expo'):
    """
    Creation of the function giving the lateral length according to its position on the parent branch

    :param pd: DataFrame - DataFrame with the laterals length law
    :param scale_x: float (0.01) - x scale by default transform x in % to real value
    :param scale_y: float (1.0) - any possible scale factor on y
    :param scale: float (1e-4) - the segment length (m)
    :param uniform: boolean or string (False) - if False use randomly an exact data point, True use a uniform distribution
            between the minimum and the maximum of the data LR_length_mm, if 'expo', use an expovariate law
    :return: a function giving the lateral length according to its position
    """
    x = pd.relative_distance_to_tip.tolist()
    y = pd.LR_length_mm.tolist()

    # size of the windows: 5%
    size = 5. * scale_x

    _length_law = histo_relative_law(x, y,
                                     size = size,
                                     scale_x = scale_x,
                                     scale_y = 1.e-3 * scale_y,
                                     scale = scale,
                                     plot = False,
                                     uniform = uniform)
    return _length_law

# to change the conductivities values by a factor to be able to do some
#    sensitivity studie
def radial(v = 92, acol = [], scale = 1):
    xr = acol[0]  # at this stage kr constant so the same x than Ka
    yr = [v * scale] * len(xr)

    return xr, yr

def axial(acol = [], scale = 1):
    x, y = acol
    y = [a * scale for a in y]
    return x, y

def my_seed():
    """ Define my own seed function to capture the seed value. """
    return int(int(_hexlify(_urandom(2500)), 16) % 100000000)

def root_creation(primary_length = 0.13, seed = None, delta = 2.0e-3, nude_length = 2.0e-2, df = None, segment_length = 1.0e-4,
                  length_data = None, branching_variability = 0.25, order_max = 4.0, order_decrease_factor = 0.7,
                  ref_radius = 7.0e-5):
    """
    creation of an mtg with properties like radius and vertex length set.

    The MTG is either generated or created from a data.
    The radius and vertex length properties are set.
    The following properties are computed: length, position, mylength, surface, volume, total length,
        primary root length, nb of intercepts

    :param:
        primary_length: primary root length for generated mtg
        seed:  seed for generated mtg, if None randomly generated
        delta: branching delay  for generated mtg
        nude_length: length from tip without lateral for generated mtg
        df: pandas DataFrame with the architecture data to be reconstructed
        segment_length: float (1.0e-4) - MTG segment length
        length_data: string - length laws file names
        branching_variability: float (0.25) - random variability of delta
        order_max: float (4.0) - maximum lateral order
        order_decrease_factor: float (0.7) - radius decrease factor between order
        ref_radius: float (7.0e-5) - the primary root radius
        intercepts: list (None) - list of distance from base
    :return:
        g: MTG with the different properties set or computed (see comments above),
        primary_length: primary root length (output for generated mtg)
        _length: total root length
        surface: total root surface
        intercepts: nb of intercepts at a given distance from base
        _seed: the seed used in the generator
    """
    if df is not None:
        g = mtg_from_aqua_data(df, segment_length)
        _seed = None
    else:
        # if no seed just create one
        if seed is None:
            _seed = my_seed()
        else:
            _seed = seed

        g = generate_g(_seed, length_data,
                       branching_variability, delta,
                       nude_length, primary_length, segment_length,
                       order_max)

    # compute radius property on MTG
    g = radius.ordered_radius(g, ref_radius, order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # Calculation of the distance from base of each vertex, used for cut and flow
    _mylength = {}
    for v in traversal.pre_order2(g, 1):
        pid = g.parent(v)
        _mylength[v] = _mylength[pid] + segment_length if pid else segment_length
    g.properties()['mylength'] = _mylength

    # _length is the total length of the RSA (sum of the length of all the segments)
    _length = g.nb_vertices(scale = 1) * segment_length
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    if df is not None:
        v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))
        primary_length = g.property('position')[v_base]

    return g, primary_length, _length, surface, _seed
