###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to calcul local relative fluxes on a given architecture and on a
#   single primary root with the same hydraulic characteristics
###############################################################################

######
# Imports

# VERSION = 2

import argparse
import sys

from openalea.mtg.algo import axis

from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters

from shared_functions import *

results = {}

################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################

parameter = Parameters()
# parameter.read_file('../example/parameters.yml')

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
    g, Keq, Jv_global = hydroroot_flow(g,
                                       segment_length = parameter.archi['segment_length'],
                                       k0 = k_radial,
                                       Jv = parameter.exp['Jv'],
                                       psi_e = parameter.exp['psi_e'],
                                       psi_base = parameter.exp['psi_base'],
                                       axial_conductivity_data = Kexp_axial_data,
                                       radial_conductivity_data = k_radial_data)

    return g, Keq, Jv_global

if __name__ == '__main__':
    seg_at_position = [1, 20, 40, 65, 100, 120, 125, 130, 135, 140, 145, 150, 155]  # distance from tip

    colors = ['orange', 'cyan', 'green', 'magenta', 'blue']

    outputfilename="fig-6D-RSA.csv"
    for iloop in range(2): # 1st for the root, 2d for cylinder because max_order set to 0 at the of the 1st pass

        nb_steps = len(parameter.output['axfold']) * len(parameter.output['radfold'])
        print('Simulation runs: ', nb_steps)
        print('#############################')
        print('figure 6-D')
        print(outputfilename)
        j_relat = {}
        _columns = []
        _columns.append('ax')
        j_relat['ax'] = []
        for i in seg_at_position:
            _columns.append(str(i) + ' mm')
            j_relat[str(i) + ' mm'] = []
        _columns.append('Jv')
        j_relat['Jv'] = []

        seed =parameter.archi['seed'][0]
        primary_length = parameter.archi['primary_length'][0]
        delta = parameter.archi['branching_delay'][0]
        nude_length = parameter.archi['nude_length'][0]


        g, primary_length, _length, surface, _seed = root_builder(primary_length = primary_length, seed = seed,
            delta = delta, nude_length = nude_length, segment_length = parameter.archi['segment_length'],
            length_data = parameter.archi['length_data'],  branching_variability = parameter.archi['branching_variability'],
            order_max = parameter.archi['order_max'], order_decrease_factor = parameter.archi['order_decrease_factor'],
            ref_radius = parameter.archi['ref_radius'])

        vertices_at_length = []
        v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))
        n_max = max(axis(g,v_base))

        for l in seg_at_position:
            ## only on PR
            vids = int(n_max-l*1.0e-3/parameter.archi['segment_length'])
            vertices_at_length.append([vids])

        j1 = {}
        for axfold in parameter.output['axfold']:
            for radfold in parameter.output['radfold']:
                avg_fold = axfold # the factor on winch the relative j is calculated
                other_fold = radfold # the other
                if avg_fold == 1: j1[other_fold] = []

                g, Keq, Jv = hydro_calculation(g, axfold = axfold, radfold = radfold)

                if avg_fold == 1:
                    g.add_property('j_relat')
                    g_1 = g.copy()
                else:
                    for v in g:
                        if v>0: g.property('j_relat')[v] = g.property('J_out')[v]/g_1.property('J_out')[v]

                c = 0
                for l in seg_at_position:
                    c += 1
                    jtot = 0.0
                    n = len(vertices_at_length[c-1])
                    for v in vertices_at_length[c-1]:
                        # remark: when done on the PR there is only 1 vertex
                        jtot += g.property('J_out')[v]

                    if avg_fold == 1:
                        j1[other_fold].append(jtot)
                        j_relat[str(l) + ' mm'].append(l*1e-3)
                    else:
                        j_relat[str(l) + ' mm'].append(jtot/j1[other_fold][c-1])

                if avg_fold == 1:
                    j1[other_fold].append(Jv)
                    j_relat['Jv'].append(primary_length)
                else:
                    j_relat['Jv'].append(Jv/j1[other_fold][c])

                j_relat['ax'].append(axfold)
                nb_steps -= 1
                sys.stdout.write('\r')
                sys.stdout.write(str(nb_steps))
                sys.stdout.flush()

        parameter.archi['order_max'] = 0 # for the cylinder

        dj2 = pd.DataFrame(j_relat, columns = _columns)
        dj1 = dj2.transpose()
        if iloop == 0:
            ax = dj1.loc['1 mm':'155 mm',[0, 1, 5, 10, 15, 19]].plot.line(x=0, color = colors, legend = False)
        else:
            dj1.loc['1 mm':'155 mm',[0, 1, 5, 10, 15, 19]].plot.line(x=0, ax = ax, style = '--', color = colors, legend = False)
            ax.set_xlabel('Distance to tip (m)')
            ax.set_ylabel('Normelized local flow (J)')
            ax.set_title('figure 6-D')
        dj1.to_csv(outputfilename, index = False, header = False)
        outputfilename = "fig-6D-cylindric.csv"
