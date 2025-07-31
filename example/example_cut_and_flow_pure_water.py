"""
Perform direct simulations or parameters adjustment to fit data of cut and flow experiment.
Water transport only, electrical network analogy

Remark:
    - Use input data see below
    - Use mainy global variables

Usage:
    %run adjustment_K_and_k.py [-h] [-o OUTPUTFILE] [-op] inputfile
        positional arguments:
          inputfile             yaml input file
        optional arguments:
          -h, --help            show this help message and exit
          -o OUTPUTFILE, --outputfile OUTPUTFILE
                                output csv file
          -op, --optimize       optimize K and k

Inputs:
    - yaml file given in command line argument
    - data/cnf_data.csv: may be changed see begining of main, csv file containing data of cut and flow data of with
            the following columns:
            - arch: sample name, the string must be contained in the 'input_file' name given in the yaml file
            - dP_Mpa: column with the working cut and flow pressure (in relative to the base) if constant, may be empty see below
            - J0, J1, ..., Jn: columns that start with 'J' containing the flux values, 1st the for the full root, then 1st cut, 2d cut, etc.
            - lcut1, ...., lcutn: columns starting with 'lcut' containing the maximum length to the base after each cut, 1st cut, 2d cut, etc. (not the for full root)
            - dP0, dP1,.., dPn: column starting with 'dP' containing the working pressure (in relative to the base) of each steps (if not constant): full root, 1st cut, 2d cut, etc.

Remark: at this stage 2022-07-29, this script is used for arabidopsis and for experiment done at a constant working pressure
        given in the yaml file, unlike adjustment_K_k_Js_Ps.py where the script has been used with CnF experiment where
        pressure may change with cut steps

Outputs:
    - console (if verbose):
        - CnF: plant name, max length (m), k (10-8 m/s/MPa), total length (m), surface (m2), Jv (microL/s)
    - matplotlib:
        - 2 plots:
            - Jv(l) cnf): Jv exp dot, Jv sim line
            - K(x): K 1st, K adjusted (displayed if adjustment asked)

    - outputfile (csv):
        - column names: 'plant', 'cut length (m)', 'primary_length (m)', 'k (10-8 m/s/MPa)', '_length (m)',
                        'surface (m2)', 'Jv (uL/s)', 'Jexp (uL/s)'
        - if Flag_Optim add the following: 'x', 'K 1st', 'K optimized'
                                            the initial and adjusted K(x)

"""

import argparse
import time
import pandas as pd
import matplotlib.pyplot as plt

from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.solver_wrapper import pure_hydraulic_model

start_time = time.time()

parser = argparse.ArgumentParser(description='run direct simulation of pure hydraulic HydroRoot, or adjust parameters on Cut and flow data.')
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("-o", "--outputfile", default = 'out.csv', help="output csv file")
parser.add_argument("-op", "--optimize", help="parameters to optimize, space separated strings, K k, "
                                              "(default: K k)", nargs='*')
parser.add_argument("-v", "--verbose", help="display convergence", action="store_true")
args = parser.parse_args()
filename = args.inputfile
output = args.outputfile
Flag_Optim = args.optimize
Flag_verbose = args.verbose
if Flag_verbose is None: Flag_verbose = False

parameter = Parameters()
parameter.read_file(filename)

### Cut and Flow DATA 
fn = 'data/maize_cnf_data.csv'
# fn = 'data/tomato_cnf_data.csv'
# fn = 'data/arabido_cnf_data.csv'
df_exp = pd.read_csv(fn, sep = ';', keep_default_na = True)
if df_exp.shape[1] == 1:
    df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)
if df_exp.shape[1] == 1:
    df_exp = pd.read_csv(fn, sep = '\t', keep_default_na = True)

### Uncomment the line below if you want to do a run using psi_base en psi_ext given in the yaml parameter file
# df_exp = None

### The run
dresults, g = pure_hydraulic_model(parameter = parameter, df_exp = df_exp, Data_to_Optim = Flag_Optim, output = output,
                                   Flag_verbose = Flag_verbose, Flag_radius = False, Flag_Constraint = False,
                                   dK_constraint = 0.0)

### Display the plot J vs Lcut
ax = dresults.plot.scatter('cut length (m)', 'Jexp (uL/s)', c = 'black')
dresults.plot.line('cut length (m)', 'Jv (uL/s)', c = 'purple', ax = ax)

### Plot K vs x and comparing radial k between 1st guess and optim value
ax_K = dresults.plot.line('x', 'K 1st', c = 'black')
dresults.plot.line('x', 'K optimized', c = 'purple', ax = ax_K)

d = pd.DataFrame({'lab':['k', 'k adjusted'], 'val':[parameter.hydro['k0'], dresults['k (10-9 m/s/MPa)'][0]]})
d.plot.bar(x='lab', y='val', rot=0)

plt.show(block=False)

