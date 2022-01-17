###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to determinate the radial conductivity k on the set of
#   architecture given in the yaml file
###############################################################################

######
# Imports

# VERSION = 2

import glob
import argparse
import sys

from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters

from shared_functions import *

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

    Flag_Optim = True

    filename = []
    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))

    fn = 'data/arabido_data.csv'
    df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)

    # predict the number of simulation run
    nb_steps = len(filename)
    print('Simulation runs: ', nb_steps)
    print('#############################')

    columns = ['plant', 'primary_length (m)', 'k (10-8 m/s/MPa)', 'total length (m)', 'surface (m2)', 'Jv (uL/s)']

    axfold = parameter.output['axfold'][0]
    radfold = parameter.output['radfold'][0]

    results = {}
    for key in columns:
        results[key] = []

    for f in filename:
        df = read_archi_data(f) if parameter.archi['read_architecture'] else None

        index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")
        sys.stdout.write('\r')
        sys.stdout.write(str(index))
        sys.stdout.flush()

        _archi_name = str(index).lower()
        j = None
        for key in df_exp['arch']:
            if key in _archi_name:
                j = df_exp[df_exp.arch == key].iloc[0].Jbase
                break
        if j is None:
            print('************ ', _archi_name, 'not in the data base')
            break
        parameter.exp['Jv'] = j

        g, primary_length, _length, surface, seed = root_creation(df = df, segment_length = parameter.archi['segment_length'],
            order_decrease_factor = parameter.archi['order_decrease_factor'], ref_radius = parameter.archi['ref_radius'])

        axial_data = parameter.hydro['axial_conductance_data']
        k0 = parameter.hydro['k0']
        Jv = 0.0

        g, Keq, Jv = hydro_calculation(g, axial_data = axial_data, k_radial = k0, axfold = axfold, radfold = radfold)

        if Flag_Optim:
            k0_old = k0
            F_old = (Jv - parameter.exp['Jv'])**2.0
            k0 *= 0.9
            eps = 1e-9
            F = 1.
            # Newton-Raphson loop to get k0
            while (F > eps):
                g, Keq, Jv = hydro_calculation(g, axial_data = axial_data, k_radial = k0, axfold = axfold, radfold = radfold)

                F = (Jv - parameter.exp['Jv']) ** 2.0

                if abs(F) > eps:
                    dfdk0 = (F - F_old) / (k0 - k0_old)

                    k0_old = k0

                    k0 = k0_old - F / dfdk0
                    while k0 < 1.0e-3:
                        k0 = 0.5 * k0_old

                    F_old = F

        results['plant'].append(index)
        results['primary_length (m)'].append(primary_length)
        results['k (10-8 m/s/MPa)'].append(k0*0.1) #uL/s/MPa/m2 -> 10-8 m/s/MPa
        results['total length (m)'].append(_length)
        results['surface (m2)'].append(surface)
        results['Jv (uL/s)'].append(Jv)

    dr = pd.DataFrame(results, columns = columns)
    print(dr.loc[:, ['plant', 'total length (m)', 'surface (m2)', 'k (10-8 m/s/MPa)']])
    if output is not None:
        dr.to_csv(output,  index = False)


