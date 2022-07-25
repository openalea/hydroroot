###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to perform some sensibility analysis on the axial and radial conductivity using
#   the architecture parameters in the yaml file
#   The factor axfold on axial data and radfold on radial k given in the parameter yaml file are used
###############################################################################

######
# Imports

# VERSION = 2

import numpy as np
import argparse
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters

from shared_functions import *

results = {}
Jv_global = 1.0

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
parameter.read_file(filename)

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None, cut_and_flow = False):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    Kexp_axial_data = axial(axial_data, axfold)
    k_radial_data = radial(k_radial, axial_data, radfold)

    # compute local jv and psi, global Jv, Keq
    g, Keq, Jv_global = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                       k0 = k_radial,
                                       Jv = parameter.exp['Jv'],
                                       psi_e = parameter.exp['psi_e'],
                                       psi_base = parameter.exp['psi_base'],
                                       axial_conductivity_data = Kexp_axial_data,
                                       radial_conductivity_data = k_radial_data)

    return g, Keq, Jv_global

if __name__ == '__main__':

    # predict the number of simulation run
    nb_steps = len(parameter.output['axfold'])*len(parameter.output['radfold'])
    print('Simulation runs: ', nb_steps)
    print('#############################')


    columns = ['seed', 'primary_length (m)', 'k (10-8 m/s/MPa)', 'ax', 'length (m)', 'surface (m2)', 'Jv (uL/s)']
    for key in columns:
        results[key] = []

    seed = parameter.archi['seed'][0]
    primary_length = parameter.archi['primary_length'][0]
    delta = parameter.archi['branching_delay'][0]
    nude_length = parameter.archi['nude_length'][0]

    g, primary_length, _length, surface, _seed = root_builder(primary_length = primary_length, seed = seed,
        delta = delta, nude_length = nude_length, df = None, segment_length = parameter.archi['segment_length'],
        length_data = parameter.archi['length_data'],  branching_variability = parameter.archi['branching_variability'],
        order_max = parameter.archi['order_max'], order_decrease_factor = parameter.archi['order_decrease_factor'],
        ref_radius = parameter.archi['ref_radius'])

    # sensibility analyse using multiplying factor on K and k
    for axfold in parameter.output['axfold']:
        for radfold in parameter.output['radfold']:
            g, Keq, Jv = hydro_calculation(g, axfold = axfold, radfold = radfold)
            results['Jv (uL/s)'].append(Jv)
            results['seed'].append(str(seed))
            results['primary_length (m)'].append(primary_length)
            results['k (10-8 m/s/MPa)'].append(parameter.hydro['k0'] * radfold * 0.1)  # uL/s/MPa/m2 -> 10-8 m/s/MPa
            results['length (m)'].append(_length)
            results['surface (m2)'].append(surface)
            results['ax'].append(axfold)


            # print nb_steps
            sys.stdout.write('\r')
            sys.stdout.write(str(nb_steps))
            sys.stdout.flush()
            nb_steps -= 1

    dresults = pd.DataFrame(results, columns = columns)
    if output is not None: dresults.to_csv(output, index = False)

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_title('figure 6-A')
    y = np.array(dresults['ax'])
    x = np.array(dresults['k (10-8 m/s/MPa)'])
    z = np.array(dresults['Jv (uL/s)'])
    scatter = ax.scatter(x,y,z)
    ax.set_xlabel(columns[2])
    ax.set_ylabel(columns[3])
    ax.set_zlabel(columns[6])

    ax2 = dresults.plot.scatter('ax', 'Jv (uL/s)', c='black')
    ax2.set_ylim([0, 0.03])
    ax2.set_title('supplemental figure 4-A')
    dresults[dresults['k (10-8 m/s/MPa)'] == 32.76558].plot.line('ax', 'Jv (uL/s)', ax = ax2, c = 'orange')
    ax2.get_legend().remove()
    ax3 = dresults.plot.scatter('k (10-8 m/s/MPa)', 'Jv (uL/s)', c='black')
    dresults[dresults['ax'] == 1.0].plot.line('k (10-8 m/s/MPa)', 'Jv (uL/s)', ax = ax3, c = 'orange')
    ax3.get_legend().remove()
    ax3.set_ylim([0, 0.03])
    ax3.set_title('supplemental figure 4-B')
