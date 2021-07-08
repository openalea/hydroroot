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

import pandas as pd
import time
import argparse

#import matplotlib.pyplot as plt
#import matplotlib as mpl
#from matplotlib.colors import Normalize

from openalea.mtg import traversal

from hydroroot import radius, markov
from hydroroot.law import histo_relative_law, reference_relative_law
from hydroroot.analysis import intercept
from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters  # import work in progress for reading init file


ONE_LAW = False
EXPOVARIATE = True

results = {}

################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################

parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("-o", "--outputfile", help="output csv file")
args = parser.parse_args()
filename = args.inputfile
output = args.outputfile
if output is None: output = "generated_architectures.csv"
parameter.read_file(filename)


def radial(v = 92, acol = [], scale = 1):
    xr = acol[0]  # at this stage kr constant so the same x than Ka
    yr = [v * scale] * len(xr)
    return xr, yr

def axial(acol = [], scale = 1):
    x, y = acol
    y = [a * scale for a in y]
    return x, y

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

    # 'expo', True, False
    _uniform = 'expo'

    # Just use the same order1 law
    law_order1 = length_law(length_data[0], scale_x = primary_length / 100., scale = segment_length, uniform = _uniform)
    if EXPOVARIATE:
        law_order2 = length_law(length_data[1], scale_x = length_max_secondary / 100., scale = segment_length,
                                uniform = _uniform)
    else:
        law_order2 = length_law(length_data[1], scale_x = length_max_secondary / 100., scale = segment_length)

    g = markov.markov_binary_tree(
        nb_vertices = nb_vertices,
        branching_variability = branching_variability,
        branching_delay = branching_delay,
        length_law = [law_order1, law_order2] if not ONE_LAW else law_order1,
        nude_tip_length = nb_nude_vertices,
        order_max = order_max,
        seed = seed)
    return g

###############################################################################
# data
###############################################################################
def length_law(pd, scale_x = 1 / 100., scale_y = 1., scale = 1e-4, uniform = True):
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

def ref_length_law(pd, scale_x = 1 / 100., scale_y = 1., scale = 1e-4, uniform = True):
    """
    scale
    """
    x = pd.relative_distance_to_tip.tolist()
    y = pd.LR_length_mm.tolist()

    # size of the windows: 5%
    size = 5. * scale_x

    _length_law = reference_relative_law(x, y,
                                         size = size,
                                         scale_x = scale_x,
                                         scale_y = 1.e-3 * scale_y)
    return _length_law

###############################################################################
# Main simulation function
###############################################################################


def my_seed():
    """ Define my own seed function to capture the seed value. """
    return int(long(_hexlify(_urandom(2500)), 16) % 100000000)

def root_creation(primary_length, seed = None, delta = 2.0e-3, nude_length = 2.0e-2, df = None):
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
    :return:
        g: MTG with the different properties set or computed (see comments above),
        primary_length: primary root length (output for generated mtg)
        _length: total root length
        surface: total root surface
        intercepts: nb of intercepts at a given distance from base
        _seed: the seed used in the generator
    """

    length_data = parameter.archi['length_data']
    g = generate_g(_seed, length_data,
                   parameter.archi['branching_variability'], delta,
                   nude_length, primary_length, parameter.archi['segment_length'],
                   parameter.archi['order_max'])

    # compute radius property on MTG
    g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])
    # compute length property and parametrisation
    g = radius.compute_length(g, parameter.archi['segment_length'])
    g = radius.compute_relative_position(g)


    # Calculation of the distance from base of each vertex, used for cut and flow
    # Remark: this calculation is done in flux.segments_at_length; analysis.nb_roots but there is a concern with the
    # parameter dl which should be equal to vertex length but which is not pass
    _mylength = {}
    for v in traversal.pre_order2(g, 1):
        pid = g.parent(v)
        _mylength[v] = _mylength[pid] + parameter.archi['segment_length'] if pid else parameter.archi['segment_length']
    g.properties()['mylength'] = _mylength

    # _length is the total length of the RSA (sum of the length of all the segments)
    _length = g.nb_vertices(scale = 1) * parameter.archi['segment_length']
    g, surface = radius.compute_surface(g)
    # g, volume = radius.compute_volume(g)


    # Compute the intercepts
    intercepts = intercept(g, sorted(parameter.output['intercepts']))


    return g, primary_length, surface, _length, intercepts, _seed

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    Kexp_axial_data = axial(axial_data, axfold)
    k_radial_data = radial(k_radial, axial_data, radfold)

    # compute local jv and psi, global Jv, Keq
    g, Keq, Jv = hydroroot_flow(g,
                                       segment_length = parameter.archi['segment_length'],
                                       k0 = k_radial,
                                       Jv = parameter.exp['Jv'],
                                       psi_e = parameter.exp['psi_e'],
                                       psi_base = parameter.exp['psi_base'],
                                       axial_conductivity_data = Kexp_axial_data,
                                       radial_conductivity_data = k_radial_data)

    return g, Keq, Jv


if __name__ == '__main__':
    count = 0

    start_time = time.time()

    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']
    run_nb = parameter.output['run_nb']
    # predict the number of simulation run
    nb_steps = run_nb * len(parameter.archi['branching_delay']) * len(parameter.output['axfold']) * \
               len(parameter.output['radfold']) * len(parameter.archi['nude_length']) * \
               len(parameter.archi['primary_length']) * len(filename)
    print 'Simulation runs: ', nb_steps
    print '#############################'


    columns = ['seed', 'primary_length (m)', 'k (10-8 m/s/MPa)', '_length (m)', 'surface (m2)', 'Jv (uL/s)']
    results = {}
    for key in columns:
        results[key] = []

    print 'seed', 'primary_length (m)', '_length (m)', 'surface (m2)', 'Jv (uL/s)'
    for i in range(run_nb):
        print run_nb - i, count
        for primary_length in parameter.archi['primary_length']:
            for delta in parameter.archi['branching_delay']:
                for nude_length in parameter.archi['nude_length']:
                    _seed = my_seed()
                    count += 1
                    g, primary_length, _length, surface, intercepts, _seed = root_creation(
                        primary_length = primary_length,
                        seed = _seed,
                        delta = delta,
                        nude_length = nude_length, df = None)
                    g, Keq, Jv = hydro_calculation(g)
                    
                    results['seed'].append(_seed)
                    results['primary_length (m)'].append(primary_length)
                    results['k (10-8 m/s/MPa)'].append(parameter.hydro['k0']*0.1) #uL/s/MPa/m2 -> 10-8 m/s/MPa
                    results['_length (m)'].append(_length)
                    results['surface (m2)'].append(surface)
                    results['Jv (uL/s)'].append(Jv)
                    print _seed, primary_length, _length, surface, Jv
    df = pd.DataFrame(results, columns = columns)
    df.to_csv(output, index = False)
    
    print time.time() - start_time
