###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to perform some sensibility analysis on the axial and radial conductivity using
#   a set of known architectures generated-roots-20-10-07.csv by there seed, primary length, internode
#   length and nude length.
#   The factor axfold on axial data and radfold on radial k given in the parameter yaml file are used
###############################################################################

######
# Imports

# VERSION = 2

import sys
import argparse
import time
import pandas as pd

from hydroroot.main import hydroroot_flow, root_builder
from hydroroot.init_parameter import Parameters
from hydroroot.conductance import axial, radial

results = {}


start_time = time.time()


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
    g, Keq, Jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                       k0 = k_radial,
                                       Jv = parameter.exp['Jv'],
                                       psi_e = parameter.exp['psi_e'],
                                       psi_base = parameter.exp['psi_base'],
                                       axial_conductivity_data = Kexp_axial_data,
                                       radial_conductivity_data = k_radial_data)

    return g, Keq, Jv

if __name__ == '__main__':

    k0 = parameter.hydro['k0']
    # dseeds = pd.read_csv('data/generated-roots-20-10-07.csv') # full set
    dseeds = pd.read_csv('data/short-generated-roots-20-10-07.csv') # short subset
    # dseeds = pd.read_csv('data/test.csv')

    # predict the number of simulation run
    nb_steps = len(dseeds) * len(parameter.output['radfold']) * len(parameter.output['axfold'])
    print('Simulation runs: ', nb_steps)
    print('#############################')

    columns = ['seed', 'primary_length (m)', 'k (10-8 m/s/MPa)', 'ax', 'length (m)', 'surface (m2)', 'Jv (uL/s)',
               'internode (m)', 'nude length (m)']
    for key in columns:
        results[key] = []

    ### variation of axfold => axial data
    count = 0
    nb_steps = len(dseeds) * len(parameter.output['axfold'])
    for id in dseeds.index:
        seed = dseeds.seed[id]
        primary_length = dseeds.primary_length[id]
        delta = dseeds.delta[id]
        nude_length = dseeds.nude_length[id]

        g, primary_length, _length, surface, _seed = root_builder(primary_length = primary_length, seed = seed,
            delta = delta, nude_length = nude_length, df = None, segment_length = parameter.archi['segment_length'],
            length_data = parameter.archi['length_data'],  branching_variability = parameter.archi['branching_variability'],
            order_max = parameter.archi['order_max'], order_decrease_factor = parameter.archi['order_decrease_factor'],
            ref_radius = parameter.archi['ref_radius'])
            

        for axfold in parameter.output['axfold']:
            g, Keq, Jv = hydro_calculation(g, axfold = axfold, k_radial = k0)
            results['Jv (uL/s)'].append(Jv)
            results['seed'].append(str(seed))
            results['primary_length (m)'].append(primary_length)
            results['k (10-8 m/s/MPa)'].append(k0 * 0.1)  # uL/s/MPa/m2 -> 10-8 m/s/MPa
            results['length (m)'].append(_length)
            results['surface (m2)'].append(surface)
            results['ax'].append(axfold)
            results['internode (m)'].append(dseeds.delta[id])
            results['nude length (m)'].append(dseeds.nude_length[id])

            
            count += 1
            progress = float(count) / nb_steps
            # print progress*100, ' %'
            sys.stdout.write('\r')
            # sys.stdout.write("[%-100s] %d%%" % ('=' * int(progress*100), int(progress*100)))
            sys.stdout.write("fig-6-C: runs done " + '{:0.4}'.format(progress*100) + ' %')
            sys.stdout.flush()

    dresults = pd.DataFrame(results, columns = columns)
    ax = dresults.plot.scatter('ax', 'Jv (uL/s)', c='black')
    ax.set_ylim([0, 0.05])
    ax.set_title('figure 6-C')

    ### variation of radfold => radial data
    for key in columns:
        results[key] = []
    count = 0
    axfold = 1.0
    nb_steps = len(dseeds) * len(parameter.output['radfold'])
    for id in dseeds.index:
            seed = dseeds.seed[id]
            primary_length = dseeds.primary_length[id]
            delta = dseeds.delta[id]
            nude_length = dseeds.nude_length[id]

            g, primary_length, _length, surface, _seed = root_builder(primary_length = primary_length, seed = seed,
                delta = delta, nude_length = nude_length, df = None, segment_length = parameter.archi['segment_length'],
                length_data = parameter.archi['length_data'],  branching_variability = parameter.archi['branching_variability'],
                order_max = parameter.archi['order_max'], order_decrease_factor = parameter.archi['order_decrease_factor'],
                ref_radius = parameter.archi['ref_radius'])

            for radfold in parameter.output['radfold']:
                g, Keq, Jv = hydro_calculation(g, axfold = axfold, radfold = radfold)
                results['Jv (uL/s)'].append(Jv)
                results['seed'].append(str(seed))
                results['primary_length (m)'].append(primary_length)
                results['k (10-8 m/s/MPa)'].append(k0 * 0.1 * radfold)  # uL/s/MPa/m2 -> 10-8 m/s/MPa
                results['length (m)'].append(_length)
                results['surface (m2)'].append(surface)
                results['ax'].append(axfold)
                results['internode (m)'].append(dseeds.delta[id])
                results['nude length (m)'].append(dseeds.nude_length[id])

                count += 1
                progress = float(count) / nb_steps
                # print progress*100, ' %'
                sys.stdout.write('\r')
                # sys.stdout.write("[%-100s] %d%%" % ('=' * int(progress*100), int(progress*100)))
                sys.stdout.write("fig-6-B: runs done " + '{:0.4}'.format(progress * 100) + ' %')
                sys.stdout.flush()

    dresults2 = pd.DataFrame(results, columns = columns)
    ax2 = dresults2.plot.scatter('k (10-8 m/s/MPa)', 'Jv (uL/s)', c = 'black')
    ax2.set_ylim([0, 0.05])
    ax2.set_title('figure 6-B')

    dr = pd.merge(dresults, dresults2, how = 'outer')
    if output is not None: dr.to_csv(output, index = False)
    print('running time is ', time.time() - start_time)