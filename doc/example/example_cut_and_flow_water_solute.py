"""
Perform direct simulations or parameters adjustment to fit data  using the water and solute transport module of Hydroroot.
It simulates two set of data: Jv(P) (flux vs pressure) and cnf (cut and flow experiment).
It may adjust parameters on either Jv(P), cnf or both data.

Remark:
    - Use input data see below

Usage:
    %run adjustment_K_k_Js_Ps.py [-h] [-o OUTPUTFILE] [-op [OPTIMIZE [OPTIMIZE ...]]] [-v] [-d DATA] inputfile

    positional arguments:
      inputfile             yaml input file
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUTFILE, --outputfile OUTPUTFILE
                            output csv file (default: out.csv)
      -op [OPTIMIZE [OPTIMIZE ...]], --optimize [OPTIMIZE [OPTIMIZE ...]]
                            parameters to optimize, space separated strings, K k
                            Ps Js sigma Klr klr, (default: K k Ps Js)
      -v, --verbose         display parameter values during adjustment (default:
                            False)
      -d [DATA], --data [DATA] data to fit: cnf, JvP or all (default: all)

Inputs:
    - yaml file given in command line argument
    - data/maize_cnf_data.csv: may be changed see begining of main, csv file containing data of cut and flow data of with
            the following columns:
            - arch: sample name that must be contained in the 'input_file' of the yaml file
            - dP_Mpa: column with the working cut and flow pressure (in relative to the base) if constant, may be empty see below
            - J0, J1, ..., Jn: columns that start with 'J' containing the flux values, 1st the for the full root, then 1st cut, 2d cut, etc.
            - lcut1, ...., lcutn: columns starting with 'lcut' containing the maximum length to the base after each cut, 1st cut, 2d cut, etc. (not the for full root) 
            - dP0, dP1,.., dPn: column starting with 'dP' containing the working pressure (in relative to the base) of each steps (if not constant): full root, 1st cut, 2d cut, etc.
    - data/maize_JvP_data.csv: may be changed see begining of main, csv file containing data of Jv(P) data of with
            the following columns:
            - arch: sample name that must be contained in the 'input_file' of the yaml file
            - J0, J1, ..., Jn: columns that start with 'J' containing the flux values of each pressure steps
            - dP0, dP1,.., dPn: column starting with 'dP' containing the working pressure (in relative to the base) of each steps

Outputs:
    - console:
        - Jv(P): DP, JvP, Cmin, Cmax, Cbase
        - CnF: max length, JvP, Cmin, Cmax, Cbase
    - matplotlib:
        - 3 subplots:
            - Jv(P): Jv exp dot, Jv sim line
            - Jv(l) cnf): Jv exp dot, Jv sim line
            - K(x): K 1st dot, K adjusted line
    - outputfile (csv):
        - column names: 'max_length', 'Jexp cnf', 'Jv cnf', 'surface', 'length', 'dp', 'Jexp(P)', 'Jv(P)', 'Cbase',
                        'kpr', 'klr', 'Js', 'Ps', 'F cnf','F Lpr', 'x pr', 'K pr', 'x lr', 'K lr',
                        'x pr', 'K1st pr', 'x lr', 'K1st lr'
        i.e.: max length from the cut to the base, J cnf exp, J cnf sim, root surface, total root length, pressure (Jv(P),
        J exp Jv(P), J sim Jv(P), solute concentration at the base (Jv(P)), radial cond PR, radial cond LR, pumping rate,
        permeability, objective fct cnf, objective fct Jv(P), x PR, K PR, x LR, K LR, x PR, K1st PR, x LR, K1st LR

"""

import argparse
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.solver_wrapper import water_solute_model

results = {}
g_cut = {}
tip_id = {}

S_g = []
cut_n_flow_length = []
Jexp = []
result_cv = []

start_time = time.time()

###############################################################################
# Command line arguments
###############################################################################


parser = argparse.ArgumentParser(description='run direct simulation of water-solute HydroRoot, or adjust parameters on Jv(P) or Cut and flow or both data.')
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("-o", "--outputfile", default = 'out.csv', help="output csv file (default: out.csv)")
parser.add_argument("-op", "--optimize", help="parameters to optimize, space separated strings, K k Ps Js sigma Klr klr, "
                                              "(default: K k Ps Js)", nargs='*')
parser.add_argument("-v", "--verbose", default = False, help="display parameter values during adjustment (default: False)", action="store_true")
parser.add_argument("-d", "--data", help="data to fit: cnf, JvP or all (default: all)", default = 'all', const = 'all', nargs='?')
args = parser.parse_args()
filename = args.inputfile
output = args.outputfile
Flag_Optim = args.optimize
Flag_verbose = args.verbose
Flag_data_to_use = args.data
if Flag_data_to_use is None: Flag_data_to_use = 'all'

parameter = Parameters()
parameter.read_file(filename)

### Cut and Flow DATA
# fn = 'data/tomato_cnf_data.csv'
fn = 'data/maize_cnf_data.csv'
df_exp = pd.read_csv(fn, sep = ';', keep_default_na = True)
if df_exp.shape[1] == 1:
    df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)
if df_exp.shape[1] == 1:
    df_exp = pd.read_csv(fn, sep = '\t', keep_default_na = True)

### Jv(P) DATA
# fn = 'data/tomato_Lpr_data.csv'
fn = 'data/maize_JvP_data.csv'
df_exp2 = pd.read_csv(fn, sep = ';', keep_default_na = True)
if df_exp2.shape[1] == 1:
    df_exp2 = pd.read_csv(fn, sep = ',', keep_default_na = True)
if df_exp2.shape[1] == 1:
    df_exp2 = pd.read_csv(fn, sep = '\t', keep_default_na = True)


### Uncomment the line below if you want to do a run using psi_base en psi_ext given in the yaml parameter file and run one JvP data point
# df_exp = df_exp2 = None
dresults, g = water_solute_model(parameter = parameter, df_cnf = df_exp, df_JvP = df_exp2,
                                Data_to_Optim = Flag_Optim, Flag_verbose = Flag_verbose, optim_method = 'COBYLA',
                                Flag_debug = Flag_verbose, Flag_Constraint = True, output = output,
                                dK_constraint = -0.03, data_to_use = Flag_data_to_use)

# 4 plots in one
################
dresults = dresults.replace(r'^s*$', float('NaN'), regex = True)
plt.ion()
fig, axs = plt.subplots(2, 2)
# Jv(P) data and fit
if 'Jexp(P)' in list(dresults.columns):
    d = dresults[['dp', 'Jexp(P)', 'Jv(P)']].dropna()
    d.sort_values(['dp'], inplace=True)
    d.plot.scatter('dp', 'Jexp(P)', ax = axs[0, 0], label = 'Jexp(P)')
    d.plot.line('dp', 'Jv(P)', ax = axs[0, 0], label = 'Jv(P)')
    j = np.array(d.loc[:, ['Jv(P)', 'Jexp(P)']])
    axs[0, 0].set_ylim(j.min(),j.max())
    axs[0,0].set_xlabel('P (Mpa)')
    axs[0,0].set_ylabel('Jv (uL/s)')

#Jv CnF data and fit
if 'Jexp cnf (uL/s)' in list(dresults.columns):
    d = dresults[['max_length', 'Jexp cnf (uL/s)', 'Jv cnf (uL/s)']].dropna()
    d.plot.scatter('max_length', 'Jexp cnf (uL/s)', ax = axs[0, 1], label = 'Jexp(P) cnf')
    d.plot.line('max_length', 'Jv cnf (uL/s)', ax = axs[0, 1], label = 'Jv(P)')
    j = np.array(d.loc[:, ['Jv cnf (uL/s)', 'Jexp cnf (uL/s)']])
    axs[0, 1].set_ylim(j.min(),j.max())
    axs[0,1].set_xlabel('max length (m)')
    axs[0,1].set_ylabel('Jv (uL/s)')

#K 1st guess and optim
d = dresults[['x pr', 'K1st pr', 'K pr']].dropna()
d.plot.scatter('x pr', 'K1st pr', ax = axs[1, 0], label = 'K1st')
d.plot.line('x pr', 'K pr', ax = axs[1, 0], label = 'K adjusted')
axs[1, 0].set_ylim(min(d['K1st pr'].min(), d['K pr'].min()),
                   max(d['K1st pr'].max(), d['K pr'].max()))
axs[1,0].set_xlabel('dist. to tip (m)')
axs[1,0].set_ylabel('K (10-9 m4/(s.Mpa)')

#radial k 1st guess and optim
d = pd.DataFrame({'lab':['k', 'k adjusted'], 'val':[parameter.hydro['k0'], dresults['kpr'][0]]})
d.plot.bar(x='lab', y='val', rot=0, ax = axs[1, 1])
axs[1,1].set_ylabel('k (10-9 m/(s.MPa)')
axs[1,1].xaxis.label.set_visible(False)
axs[1,1].get_legend().remove()

fig.patch.set_facecolor('lightgrey')
fig.tight_layout()

print('running time is ', time.time() - start_time)
