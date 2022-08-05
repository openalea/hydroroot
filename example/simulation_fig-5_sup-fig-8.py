###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to calcul the outgoing flux on a set of known architectures
#   generated-roots-20-10-07.csv by there seed, primary length, internode length and nude length
#   generate Figures: J vs primary length, vs nude length and internode length
#   If the argument outpufile is not None save the result to a csv file
###############################################################################

import pandas as pd
import sys
import matplotlib.pyplot as plt
import argparse
import time

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
    # dseeds = pd.read_csv('data/generated-roots-20-10-07.csv') # Complete set
    dseeds = pd.read_csv('data/short-generated-roots-20-10-07.csv') # 10 times shortest usefull for run time reduction


    # if a seed is given in the parameters.yml file then restrict to this seed
    if parameter.archi['seed'][0] is not None:
        dseeds = dseeds[dseeds.seed == parameter.archi['seed']]

    # predict the number of simulation run
    nb_steps = len(dseeds)
    print('Simulation runs: ', nb_steps)
    print('#############################')

    for i in range(5):
        results[i] = {}
        columns = ['seed', 'primary_length (m)', 'k (10-8 m/s/MPa)', 'ax', 'length (m)', 'surface (m2)', 'Jv (uL/s)',
                   'internode (m)', 'nude length (m)']
        for key in columns:
            results[i][key] = []
    count = 0
    # sensibility analyse using multiplying factor on K and k
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
            
        i = 0
        for axfold, k0 in [(0.5, 98.6), (1.0, 327.6), (2.0, 714.3), (5.0, 714.3), (0.01, 98.6)]:
            g, Keq, Jv = hydro_calculation(g, axfold = axfold, k_radial = k0)
            results[i]['Jv (uL/s)'].append(Jv)
            results[i]['seed'].append(str(seed))
            results[i]['primary_length (m)'].append(primary_length)
            results[i]['k (10-8 m/s/MPa)'].append(k0 * 0.1)  # uL/s/MPa/m2 -> 10-8 m/s/MPa
            results[i]['length (m)'].append(_length)
            results[i]['surface (m2)'].append(surface)
            results[i]['ax'].append(axfold)
            results[i]['internode (m)'].append(dseeds.delta[id])
            results[i]['nude length (m)'].append(dseeds.nude_length[id])
            i += 1
            
        count += 1
        progress = float(count) / len(dseeds)
        # print progress*100, ' %'
        sys.stdout.write('\r')
        # sys.stdout.write("[%-100s] %d%%" % ('=' * int(progress*100), int(progress*100)))
        sys.stdout.write("runs done " + '{:0.4}'.format(progress*100) + ' %')
        sys.stdout.flush()


    dresults = pd.DataFrame(results[0], columns = columns)

    #If the argument outpufile is not None save the result to a csv file
    if output is not None: dresults.to_csv(output, index = False)

    ###### Figures creation #####

    #Figures 5
    ax = dresults.plot.scatter('surface (m2)', 'Jv (uL/s)', color='Blue', edgecolors = 'Blue')
    ax.set_ylim([0, 0.06]), ax.set_xlim([0, 20e-4])
    ax.set_title('figure 5-A')

    dresults = pd.DataFrame(results[1], columns = columns)
    # if output is not None: dresults.to_csv(output, index = False)
    dresults.plot.scatter('surface (m2)', 'Jv (uL/s)', ax = ax, color='orange', edgecolors = 'orange')
    ax2 = dresults.plot.scatter('primary_length (m)', 'Jv (uL/s)')
    ax2.set_title('figure 5-B')
    ax2.axis(xmin = 0.0, xmax = max(results[1]['primary_length (m)'])*1.1, ymin = 0.0, ymax = max(results[1]['Jv (uL/s)'])*1.1)
    ax3 = dresults.plot.scatter('internode (m)', 'Jv (uL/s)')
    ax3.axis(xmin = 0.0, xmax = max(results[1]['internode (m)'])*1.1, ymin = 0.0, ymax = max(results[1]['Jv (uL/s)'])*1.1)
    ax3.set_title('figure 5-C')
    ax4 = dresults.plot.scatter('nude length (m)', 'Jv (uL/s)')
    ax4.axis(xmin = 0.0, xmax = max(results[1]['nude length (m)'])*1.1, ymin = 0.0, ymax = max(results[1]['Jv (uL/s)'])*1.1)
    ax4.set_title('figure 5-D')

    dresults = pd.DataFrame(results[2], columns = columns)
    # if output is not None: dresults.to_csv(output, index = False)
    dresults.plot.scatter('surface (m2)', 'Jv (uL/s)', ax = ax, color='green', edgecolors = 'green')

    drealplants = pd.read_csv('data/10-arabido-plants.csv')
    drealplants.plot.scatter('surface', 'Jv', ax = ax, color='black', s=25)

    #supplemental figure 6
    fig = {}
    ax5 = {}
    ax6 = {}
    for s in ['primary_length (m)', 'internode (m)', 'nude length (m)']:
        fig[s] = plt.figure()
        ax5[s] = fig[s].add_subplot(111, label = "1")
        dresults = pd.DataFrame(results[3], columns = columns)
        dresults.plot.scatter(s, 'Jv (uL/s)', ax = ax5[s], color = 'black')
        ax5[s].set_title('supplemental figure 8')
        ax5[s].axis(xmin = 0.0, xmax = max(results[3][s])*1.1, ymin = 0.0, ymax = max(results[3]['Jv (uL/s)'])*1.1)
        xmax = ax5[s].get_xlim()[1] * 0.99
        xmin = ax5[s].get_xlim()[0] - 0.01 * ax5[s].get_xlim()[1]

        ax6[s] = fig[s].add_subplot(111, label = "2", frame_on = False)
        dresults = pd.DataFrame(results[4], columns = columns)
        dresults.plot.scatter(s, 'Jv (uL/s)', ax = ax6[s], color = 'orange', edgecolors = 'orange')
        ax6[s].set_xlim(xmin,xmax)
        ax6[s].axis(ymin = 0.0, ymax = max(results[4]['Jv (uL/s)'])*1.1)
        ax6[s].yaxis.tick_right()
        ax6[s].yaxis.set_label_position('right')
        ax6[s].get_xaxis().set_visible(False)
        ax6[s].tick_params(axis = 'y', color = "orange")

    print('running time is ', time.time() - start_time)