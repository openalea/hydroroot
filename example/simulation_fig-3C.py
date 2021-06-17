###############################################################################
#
# Authors: C. Pradal, Y. Boursiac
# Date : 14/10/2016
#
# Date: 2019-12-03
# Modified by F. Bauget to test yaml configuration file
#
# Date: 2019-12-10
# F. Bauget merging simulation.py and hydro_measures
###############################################################################

######
# Imports

# VERSION = 2

from random import _hexlify, _urandom

import numpy as np
import pandas as pd
import argparse
import tempfile, os

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize

from openalea.plantgl.all import Viewer
from IPython.display import Image, display

from hydroroot import radius, markov
from hydroroot.law import histo_relative_law
from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters  # import work in progress for reading init file
from hydroroot.display import plot as mtg_scene

################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################

parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile")
args = parser.parse_args()
filename = args.inputfile
parameter.read_file(filename)

# read architecture file
def read_archi_data(fn):
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
    scale
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
    return int(long(_hexlify(_urandom(2500)), 16) % 100000000)

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None, cut_and_flow = False):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    kexp_axial_data = axial(axial_data, axfold)
    k_radial_data = radial(k_radial, axial_data, radfold)

    # compute local jv and psi, global Jv, Keq
    g, keq, jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                   psi_e = parameter.exp['psi_e'],
                                   psi_base = parameter.exp['psi_base'],
                                   axial_conductivity_data = kexp_axial_data,
                                   radial_conductivity_data = k_radial_data)

    return g, keq, jv

def plot(g, name=None, **kwds):

    Viewer.display(mtg_scene(g, **kwds))
    if name is not None:
            Viewer.frameGL.saveImage(name)

def my_plot_with_bar(prop):
    values = np.array(g.property(prop).values())
    plot(g, prop_cmap = prop, lognorm = True)
    fig, ax = plt.subplots(figsize = (6, 1))
    fig.subplots_adjust(bottom = 0.5)
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin = values.min(), vmax = values.max())
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'horizontal')
    cb1.set_label(prop)
    fig.show()

if __name__ == '__main__':

    count = 0
    for seed in parameter.archi['seed']:
        length_data = parameter.archi['length_data']
        g = generate_g(seed, length_data,
                       parameter.archi['branching_variability'], parameter.archi['branching_delay'][count],
                       parameter.archi['nude_length'][count], parameter.archi['primary_length'][count],
                       parameter.archi['segment_length'],parameter.archi['order_max'])

        # compute radius property on MTG
        g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])

        # compute length property and parametrisation
        g = radius.compute_length(g, parameter.archi['segment_length'])
        g = radius.compute_relative_position(g)

        g, Keq, Jv = hydro_calculation(g)

        # plot(g, name = 'plot-' + str(seed) + '.png', prop_cmap = 'order')
        print 'seed: ', seed
        prop_cmap = 'order'

        plot(g, prop_cmap = prop_cmap)
        fn = tempfile.mktemp(suffix = '.png')
        Viewer.saveSnapshot(fn)
        Viewer.stop()
        img = Image(fn)
        os.unlink(fn)
        display(img)

        count += 1
