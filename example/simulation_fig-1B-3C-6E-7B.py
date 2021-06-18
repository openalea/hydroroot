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
from random import _hexlify, _urandom

import pandas as pd
import glob
import argparse
import tempfile, os

import openalea.plantgl.all as pgl
from openalea.mtg import turtle as turt
from IPython.display import Image, display

from hydroroot import radius, markov
from hydroroot.law import histo_relative_law
from hydroroot.generator.measured_root import mtg_from_aqua_data
from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters
from hydroroot.display import get_root_visitor
from hydroroot.display import plot as mtg_scene

################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################

parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("--prop", help="property to display, e.g.: order, j")
args = parser.parse_args()
filename = args.inputfile
prop = args.prop
if prop is None: prop = 'order'
parameter.read_file(filename)

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

def plot_order(g1, has_radius=True, r_base=1.e-4, r_tip=5e-5, prune=None, name=None):
    """
    Display the architecture in plantGL Viewer with roots colors according to there order
    The radius property may be changed for display purpose.
    The MTG g1 stay unmodified

    :param g1: MTG() - the architecture to display
    :param has_radius: Boolean (True) - True use the radius property values, calculate them otehrwise according to r_base and r_tip
    :param  r_base: float (1e-4) - if has_radius is False, the radius at the base of a root whatever its order (mm)
    :param  r_tip: float (5e-5) - if has_radius is False, the radius at the tip of a root whatever its order (mm)
    :param  prune: float (None) - distance from the base of the primary after which the root is not displayed
    :param name: string (None) - if not None, the name of the saved file
    :return:
    """
    g = g1.copy() # because we may change the radius if we want
    visitor = get_root_visitor(prune=prune)

    # changing radius just for display
    r_base, r_tip = float(r_base), float(r_tip)
    if not has_radius:
        radius.discont_radius(g,r_base=r_base, r_tip=r_tip)

    turtle = turt.PglTurtle()
    turtle.down(180)
    scene = turt.TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False)

    c = {}
    c[0] = [0,0,255]
    c[1] = [0,127,0]
    c[2] = [255,0,0]
    c[3] = c[2]
    for v in g.vertices_iter(scale = g.max_scale()):
        o = g.property('order')[v]
        g.property('color')[v] = c[o]

    shapes = dict( (sh.getId(),sh) for sh in scene)

    colors = g.property('color')
    for vid in colors:
        if vid in shapes:
            shapes[vid].appearance = pgl.Material(colors[vid])
    scene = pgl.Scene(shapes.values())

    pgl.Viewer.display(scene)
    if name is not None:
            pgl.Viewer.frameGL.saveImage(name)


def plot(g, name=None, **kwds):
    """
    Display the architecture in plantGL Viewer with roots colors according to the property chosen
    :param g: MTG()
    :param name: string - if not None, the name of the saved file
    :param kwds: parameters of hydroroot.display.plot()
        - has_radius: Boolean (False) - True use the radius property values, calculate them otehrwise according to r_base and r_tip
        - r_base: float (1e-4) - if has_radius is False, the radius at the base of a root whatever its order (mm)
        - r_tip: float (5e-5) - if has_radius is False, the radius at the tip of a root whatever its order (mm)
        - prune: float (None) - distance from the base of the primary after which the root is not displayed
        - prop_cmap: string ('radius') - the property name used for the color map
        - cmap: string ('jet') - the name of the matplotlib colormap to use
        - lognorm: Boolean (False) - True: log-normalised, normalised otherwise
    :return:
    """
    pgl.Viewer.display(mtg_scene(g, **kwds))
    if name is not None:
            pgl.Viewer.frameGL.saveImage(name)

if __name__ == '__main__':

    # files names of the architecture if reconstructed from a file
    # if not we just give a dummy name for the loop used to launch run
    filename = []
    if parameter.archi['read_architecture']:
        run_nb = 1
        parameter.archi['seed'] = [1] # to give something to the For-Loop below
        for f in parameter.archi['input_file']:
            filename = filename + (glob.glob(parameter.archi['input_dir'] + f))
    else:
        run_nb = parameter.output['run_nb']
        filename = ['one_run']  # just to have one run in the For-Loop below
        if parameter.archi['seed'] is None:
            s = my_seed()
            parameter.archi['seed'] = list(s)


    for seed in parameter.archi['seed']:
        for f in filename:
            if parameter.archi['read_architecture']:
                df = read_archi_data(f)
                g = mtg_from_aqua_data(df, parameter.archi['segment_length'])
            else:
                length_data = parameter.archi['length_data']
                g = generate_g(seed, length_data,
                               parameter.archi['branching_variability'], parameter.archi['branching_delay'][0],
                               parameter.archi['nude_length'][0], parameter.archi['primary_length'][0],
                               parameter.archi['segment_length'],parameter.archi['order_max'])

            # compute radius property on MTG
            g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])

            # compute length property and parametrisation
            g = radius.compute_length(g, parameter.archi['segment_length'])
            g = radius.compute_relative_position(g)

            for axfold in parameter.output['axfold']:
                g, Keq, Jv = hydro_calculation(g, axfold = axfold)

                if parameter.archi['read_architecture']:
                    index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")
                else:
                    index = seed

                #prop_cmap: property to plot e.g.: for figure 1B prop_'order', for figure 3CD 'j'
                # it could also be 'J_out' the axial flux
                print index, axfold
                # prop_cmap = 'order'

                # g has radius, here we set fictive radii just for visual comfort
                alpha = 0.2 # radius in millimeter identical for all orders
                gcopy = g.copy() # copy because we change the radius property in plot below
                if prop != 'order':
                    plot(gcopy, has_radius=False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4, prop_cmap = prop)
                else:
                    plot_order(gcopy, has_radius=False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4)
                pgl.Viewer.widgetGeometry.setSize(450, 600) # set the picture size in px
                fn = tempfile.mktemp(suffix='.png')
                pgl.Viewer.saveSnapshot(fn)
                pgl.Viewer.stop()
                img = Image(fn)
                os.unlink(fn)
                display(img)