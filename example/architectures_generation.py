###############################################################################
#
# Authors: C. Pradal, Y. Boursiac
# Date : 14/10/2016
#
# Date: 2019-12-03
# Modified by F. Bauget
#
# Date: 2025-04-16
# this is just an example of several architectures generation according to parameters_architecture_generation.yml
# that presents range of values for different parameter (primary_length, branching_delay, etc.)
# so loops are done on the different parameters and a direct simulation for each architecture is done
# in parameters_architecture_generation the parameter run_nb can also be changed to do several run with the same parameter
# set (does not mean same architecture because the seeds will be different)

# "%run architectures_generation parameters_architecture_generation.yml" in a python console
###############################################################################

import pandas as pd
import time
import argparse

from hydroroot.conductance import radial
from hydroroot.main import hydroroot_flow, root_builder
from hydroroot.init_parameter import Parameters

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

if __name__ == '__main__':
    count = 0

    start_time = time.time()

    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']
    run_nb = parameter.output['run_nb']
    # predict the number of simulation run
    nb_steps = run_nb * len(parameter.archi['branching_delay']) *  \
               len(parameter.archi['nude_length']) * \
               len(parameter.archi['primary_length'])
    print('Simulation runs: ', nb_steps)
    print('#############################')


    columns = ['seed', 'primary_length (m)', 'k (10-8 m/s/MPa)', '_length (m)', 'surface (m2)', 'Jv (uL/s)']
    results = {}
    for key in columns:
        results[key] = []

    axial_data = parameter.hydro['axial_conductance_data']
    k_radial_data = radial(parameter.hydro['k0'], axial_data)

    print('seed', 'primary_length (m)', '_length (m)', 'surface (m2)', 'Jv (uL/s)')
    for i in range(run_nb):
        if run_nb > 1: print(run_nb - i, count)
        for primary_length in parameter.archi['primary_length']:
            for delta in parameter.archi['branching_delay']:
                for nude_length in parameter.archi['nude_length']:
                    count += 1
                    g, p, _length, surface, _seed = root_builder( primary_length = primary_length,
                                                               delta = delta,
                                                               nude_length = nude_length, df = None,
                                                               segment_length = parameter.archi['segment_length'],
                                                               length_data = parameter.archi['length_data'],
                                                               order_max = parameter.archi['order_max'],
                                                               order_decrease_factor = parameter.archi['order_decrease_factor'],
                                                               ref_radius = parameter.archi['ref_radius'])

                    g, Keq, Jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                                k0 = k_radial_data,
                                                psi_e = parameter.exp['psi_e'],
                                                psi_base = parameter.exp['psi_base'],
                                                axial_conductivity_data = axial_data,
                                                radial_conductivity_data = k_radial_data)
                    
                    results['seed'].append(_seed)
                    results['primary_length (m)'].append(primary_length)
                    results['k (10-8 m/s/MPa)'].append(parameter.hydro['k0']*0.1) #uL/s/MPa/m2 -> 10-8 m/s/MPa
                    results['_length (m)'].append(_length)
                    results['surface (m2)'].append(surface)
                    results['Jv (uL/s)'].append(Jv)
                    # print(_seed, primary_length, _length, surface, Jv)
    df = pd.DataFrame(results, columns = columns)
    df.to_csv(output, index = False)
    
    print(time.time() - start_time)
