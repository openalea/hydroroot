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
    - data/maize_Lpr_data.csv: may be changed see begining of main, csv file containing data of Jv(P) data of with
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
import glob
import argparse
import copy
import math
import time
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from openalea.mtg.algo import axis
from scipy import optimize, constants
# from skopt import gp_minimize #comes from version 2.7

from hydroroot.read_file import read_archi_data
from hydroroot.main import root_builder
from hydroroot import radius, flux
from hydroroot.init_parameter import Parameters
from hydroroot.conductance import set_conductances, axial
from hydroroot.water_solute_transport import pressure_calculation, pressure_calculation_no_non_permeating_solutes, \
    init_some_MTG_properties, osmotic_p_peg


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

parameter = Parameters()

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
if Flag_Optim is None:
    Flag_Optim = [] # if -op not passed then we have an empty lits
elif len(Flag_Optim) == 0:
    Flag_Optim = ['K', 'k', 'Js', 'Ps']  # if -op is passed without string
Flag_verbose = args.verbose
Flag_data_to_use = args.data
parameter.read_file(filename)

###############################################################################
# Some Global boolean Flags: Flags I don't change often could be passed by command line
###############################################################################

Flag_radius = True  # True if radii furnished in architecture file used them, other wise use ref_radius and so on
Flag_Optim_K = ('K' in Flag_Optim)  # optimize axial conductance K
Flag_Constraint = True  # use constraint dK_constraint on 1st derivative dK/dx of K
Flag_Optim_Klr = ('Klr' in Flag_Optim)  # optimize axial conductance of laterals Klr if <> than PR
Flag_Optim_klr = ('klr' in Flag_Optim)  # optimize radial conductivity of laterals klr if <> than PR
Flag_Optim_Js = ('Js' in Flag_Optim)    # optimize pumping rate Js
Flag_Optim_Ps = ('Ps' in Flag_Optim)  # optimize permeability
Flag_Optim_k = ('k' in Flag_Optim)  # optimize radial conductivity k
Flag_Optim_sigma = ('sigma' in Flag_Optim)
Flag_w_Lpr = False  # set weight to  1.0 / len(list_DP)
Flag_w_cnf = False  # set weight to  len(cut_n_flow_length)

###############################################################################
# Functions
###############################################################################

def fun_constraint(x):
    """
    Calculation of the constraint for the optimize.minimize solver
    array of 1 column with non-negative constraints, i.e. every c[i] >= 0
    :param x: array of the parameters to optimize
    :return: numpy arrays
    """
    n = len(axial_data[1])
    c = np.ones(n - 1)
    l = axial_data[0]

    # inequality constraints >= 0 for the axial conductance: line below <=> (x[i+1] - x[i])/(l[i+1] - l[i]) >= dK_constraint
    for i in range(n - 1):
        c[i] = x[i+1] - x[i] - dK_constraint * (l[i+1] - l[i])

    # bounds as inequality constraints needed by solver 'COBYLA' but redundant for 'SLSQP' with bounds
    for i in range(n,len(x)):
        c = np.append(c, (x[i] - bnds[i][0]))
        c = np.append(c, (bnds[i][1] - x[i]))

    return c

def fun_bound_cobyla(x):
    """
    Calculation of the constraint for the optimize.minimize solver COBYLA
    bounds expressed as non-negative constraint
    array of 1 column with non-negative constraints, i.e. every c[i] >= 0
    :param x: array of the parameters to optimize
    :return: numpy arrays
    """
    c = []
    for i in range(len(x)):
        c.append(x[i] - bnds[i][0])
        c.append(bnds[i][1] - x[i])
    return np.array(c)

def fun(x):
    """
    Calculation of the objective function (Residual sum of squares) done on Jv(P) and CnF data

    :param x: array of the parameters to optimize
    :return: F (float), the objective function
    """

    # Kx pr, Kx lr, Js, Ps,k lr and k pr may be optimized depends on the Flags_Optim_
    _x = x * xini
    if iKpr > 0:
        axial_data[1] = list(_x[:iKpr + 1]) # np array : _x[:n] means the n 1st elements index [0;n-1] but _x[n] means _x index n

    if not Flag_Optim_Klr:
        _axial_lr = None
    else:
        _axial_lr = axial_lr
        _axial_lr[1] = _x[iKpr + 1:iKlr + 1]

    if iJs > 0:
        Js = _x[iJs]
    else:
        Js = J_s
    if iPs > 0:
        Ps = _x[iPs]
    else:
        Ps = P_s

    if Flag_Optim_klr:
        klr = _x[iklr]
    else:
        klr = None

    if Flag_Optim_k:
        kpr = _x[ikpr]
    else:
        kpr = k[0]
    if Flag_Optim_sigma:
        sigma = _x[isig]
    else:
        sigma = Sigma
    # set new K and k in the MTG
    g = set_K_and_k(g_cut, axial_data, kpr, axial_lr = _axial_lr, k_lr = klr, nr = Nb_of_roots,
                        nl = len(cut_n_flow_length))
    # run simulation Jv(P) and CnF
    JvP, F_JvP, C = Jv_P_calculation(g, sigma, Js, Ps)
    JvCnf, F_JvCnF, C = Jv_cnf_calculation(g, sigma, Js, Ps)

    F = F_JvP + F_JvCnF

    ## with COBYLA constraints are not always respected, here a hugly way to force them by increasing F when they are not respected
    # n = len(axial_data[1])
    # c = np.ones(n - 1)
    # l = axial_data[0]
    # for i in range(n - 1):
    #     if _x[i+1] - _x[i] - dK_constraint * (l[i+1] - l[i]) < 0.0:
    #         F += 1.0

    if Flag_verbose: print('{:0.3e}'.format(F), ' '.join('{:0.3e}'.format(i) for i in _x))

    result_cv.append([F,kpr,klr,Ps,Js] + axial_data[1])
    return F

def fun_JvP_only(x):
    """
    Calculation of the objective function (Residual sum of squares) done on Jv(P) data

    :param x: array of the parameters to optimize
    :return: F (float), the objective function
    """

    # Kx pr, Kx lr, Js, Ps,k lr and k pr may be optimized depends on the Flags_Optim_
    _x = x * xini
    if iKpr > 0:
        axial_data[1] = list(_x[:iKpr + 1]) # np array : _x[:n] means the n 1st elements index [0;n-1] but _x[n] means _x index n

    if not Flag_Optim_Klr:
        _axial_lr = None
    else:
        _axial_lr = axial_lr
        _axial_lr[1] = _x[iKpr + 1:iKlr + 1]

    if iJs > 0:
        Js = _x[iJs]
    else:
        Js = J_s
    if iPs > 0:
        Ps = _x[iPs]
    else:
        Ps = P_s

    if Flag_Optim_klr:
        klr = _x[iklr]
    else:
        klr = None

    if Flag_Optim_k:
        kpr = _x[ikpr]
    else:
        kpr = k[0]

    if Flag_Optim_sigma:
        sigma = _x[isig]
    else:
        sigma = Sigma
    # set new K and k in the MTG
    g = set_K_and_k(g_cut, axial_data, kpr, axial_lr = _axial_lr, k_lr = klr, nr = Nb_of_roots, nl = len(cut_n_flow_length))
    # run simulation Jv(P)
    JvP, F, C = Jv_P_calculation(g, sigma, Js, Ps)

    if Flag_verbose: print('{:0.3e}'.format(F), ' '.join('{:0.3e}'.format(i) for i in _x))

    result_cv.append([F,kpr,klr,Ps,Js] + axial_data[1])
    return F

def fun_cnf_only(x):
    """
    Calculation of the objective function (Residual sum of squares) done on CnF data

    :param x: array of the parameters to optimize
    :return: F (float), the objective function
    """

    # Kx pr, Kx lr, Js, Ps,k lr and k pr may be optimized depends on the Flags_Optim_
    _x = x * xini
    if iKpr > 0: axial_data[1] = list(_x[:iKpr + 1])  # np array : _x[:n] means the n 1st elements index [0;n-1] but _x[n] means _x index n

    if not Flag_Optim_Klr:
        _axial_lr = None
    else:
        _axial_lr = axial_lr
        _axial_lr[1] = _x[iKpr + 1:iKlr + 1]

    if iJs > 0:
        Js = _x[iJs]
    else:
        Js = J_s
    if iPs > 0:
        Ps = _x[iPs]
    else:
        Ps = P_s

    if Flag_Optim_klr:
        klr = _x[iklr]
    else:
        klr = None

    if Flag_Optim_k:
        kpr = _x[ikpr]
    else:
        kpr = k[0]

    if Flag_Optim_sigma:
        sigma = _x[isig]
    else:
        sigma = Sigma
    # set new K and k in the MTG
    g = set_K_and_k(g_cut, axial_data, kpr, axial_lr = _axial_lr, k_lr = klr, nr = Nb_of_roots, nl = len(cut_n_flow_length))
    # run simulation Jv(P)
    JvCnf, F, C = Jv_cnf_calculation(g, sigma, Js, Ps)

    if Flag_verbose: print('{:0.3e}'.format(F), ' '.join('{:0.3e}'.format(i) for i in _x))

    result_cv.append([F, kpr, klr, Ps, Js] + axial_data[1])
    return F

def set_K_and_k(g, axial_pr, k_pr, axial_lr = None, k_lr = None, nr = 1, nl = 0):

    """
    set the axial conductance and the radial conductivity of the uncut root and the different cut roots
    The vertices of the cut roots from the entire root may have changed therefore the vertices where the cuts are made
    must be set with the correct K and k see the code

    :param g: (dict) - dictionnary of MTG corresponding to the entire root and the cuts
    :param axial_pr: (list) - the axial conductance, list of 2 lists of floats
    :param k_pr:  (float) - the radial conductivity
    :param axial_lr: (list) - if not None the axial conductance of the laterals, list of 2 lists of floats
    :param k_lr: (float) - if not None the radial conductivity  of the laterals
    :param nr: (int) - number of root, because with seminals the measurements may have been done with several roots
    :param nl: (int) - number of cuts
    :return: g
    """

    # if axial_lr is None: axial_lr = axial_pr
    # if k_lr is None: k_lr = k_pr

    for ig in range(nr):
        g[0, ig] = set_conductances(g[0, ig], axial_pr = axial_pr, k0_pr = k_pr, axial_lr = axial_lr, k0_lr = k_lr)
        # set the different cut roots
        for ic in range(1, nl + 1):
            for v in g[ic, ig].vertices_iter(g[0, ig].max_scale()):
                vid = g[ic, ig].property('original_vid')[v]
                g[ic, ig].property('K')[v] = g[0, ig].property('K')[vid]
                g[ic, ig].property('k')[v] = g[0, ig].property('k')[vid]
                g[ic, ig].property('K_exp')[v] = g[0, ig].property('K_exp')[vid] # needed in pressure_calculation if Cpeg because calculation of K
                g[ic, ig].property('k0')[v] = g[0, ig].property('k0')[vid]
                # difference with the resistance network we do not set k = K the real boundary condition is used
                # with the help of label 'cut' see pressure_calculation
    return g

def Jv_P_calculation(g, sigma, Js, Ps):
    """
    Perform the calculation of the data Jv(P), i.e. for different pressure difference list_DP
    Most of the variables are global variables, only the variables that change at each calculation g (MTG) or
    that are able to be optimized (sigma, Js, Ps) are passed in arguments
    The change of K and k have been taken into account in function set_K_and_k

    Return
    JvP a dictionnary of outgoing  sap flux at each pressure step and for each roots (for the case they are several seminals)
    C  a dictionnary of the concentration map at each pressure step and for each roots (for the case they are several seminals)
    F the objective fonction

    :param g: (dict) - dictionnary of MTG corresponding to the entire root and the cuts
    :param sigma: (float) - the reflection coefficient
    :param Js: (float) - the pumping rate
    :param Ps: (float) - the permeability
    :return: JvP (dict), C {dict}, F (float)
    """
    JvP = {}
    C = {}
    F = 0.0
    data=None
    row=None
    col=None
    C_base = C_base_ini
    for idP in range(len(list_DP)):
        Jv = 0.0
        JCbase = 0.0
        for ig in range(Nb_of_roots):
            g[0, ig] = flux.flux(g[0, ig], psi_e = psi_base + list_DP[idP], psi_base = psi_base, invert_model = True)
            g[0, ig] = init_some_MTG_properties(g[0, ig], tau = Js, Cini = Cini, t = 1, Ps = Ps)
            nb_v = g[0, ig].nb_vertices()
            Fdx = 1.0
            Fdx_old = 1.
            Jv_old = 1.
            # Newton-Raphson schemes: in pressure_calculation_no_non_permeating_solutes calculation of dx,
            # array with dP and dC variation of the variables between two Newton step. Then the Newton scheme stops when
            # Fdx > eps see below
            while Fdx > eps:
                g[0, ig], dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g[0, ig], sigma = sigma, Ce = Ce,
                                                                   Pe = parameter.exp['psi_e'], Pbase = parameter.exp['psi_base'],
                                                                    Cse = Cse, dP = list_DP[idP], C_base = None)
                Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
                JvP[idP, ig] = g[0, ig].property('J_out')[1]
                if abs(JvP[idP, ig] - Jv_old) < 1.0e-4:
                    break
                if abs(Fdx - Fdx_old) < eps:
                    break
                Fdx_old = Fdx
                Jv_old = JvP[idP, ig]

            Jv += JvP[idP, ig]

            C[idP, ig] = copy.deepcopy(g[0, ig].property('C'))
            if Jv >= 0.0: JCbase += JvP[idP, ig] * C[idP, ig][1]
        F += w_Lpr * (Jv - list_Jext[idP])**2.0 #/ S_g[0]**2.0 #/ float(len(list_DP)) #/ list_Jext[idP]**2.0
        if JCbase > 0.0: C_base = JCbase / Jv
    return JvP, F, C

def Jv_cnf_calculation(g, sigma, Js, Ps):
    """
    Perform the calculation of the data CnF, i.e. for different cut length cut_n_flow_length
    Most of the variables are global variables, only the variables that change at each calculation g (MTG) or
    that are able to be optimizes (sigma, Js, Ps) are passed in arguments

    Return
    JvCnf a dictionnary of outgoing  sap flux at each cut step and for each roots (for the case they are several seminals)
    C  a dictionnary of the concentration map at each cut step and for each roots (for the case they are several seminals)
    F the objective fonction

    :param g: (dict) - dictionnary of MTG corresponding to the entire root and the cuts
    :param sigma: (float) - the reflection coefficient
    :param Js: (float) - the pumping rate
    :param Ps: (float) - the permeability
    :return: JvCnf (dict), C {dict}, F (float)
    """
    ic = 0
    JvCnf = {}
    C = {}
    F = 0.0
    Jv = 0.0
    C_base = C_base_ini
    JCbase = 0.0
    data = row = col = None
    for ig in range(Nb_of_roots):
        g[ic, ig] = flux.flux(g[ic, ig], psi_e = psi_base + DP_cnf[ic], psi_base = psi_base, invert_model = True)
        g[ic, ig] = init_some_MTG_properties(g[ic, ig], tau = Js, Cini = Cini, t = 1, Ps = Ps)
        nb_v = g[ic, ig].nb_vertices()
        Fdx = 1.0
        Fdx_old = 1.
        Jv_old = 1.
        # Newton-Raphson schemes: in pressure_calculation_no_non_permeating_solutes calculation of dx,
        # array with dP and dC variation of the variables between two Newton step. Then the Newton scheme stops when
        # Fdx > eps see below
        while Fdx > eps:
            # use pressure_calculation_no_non_permeating_solutes because the root is uncut so no PEG enter the root
            g[ic, ig], dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g[ic, ig], sigma = sigma,
                                                                                           Ce = Ce,
                                                                                           Pe = parameter.exp['psi_e'],
                                                                                           Pbase = parameter.exp[
                                                                                               'psi_base'],
                                                                                           Cse = Cse,
                                                                                           dP = DP_cnf[ic])
            Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
            JvCnf[ic, ig] = g[ic, ig].property('J_out')[1]
            if abs(JvCnf[ic, ig] - Jv_old) < 1.0e-4:
                break
            if abs(Fdx - Fdx_old) < eps:
                break
            Fdx_old = Fdx
            Jv_old = JvCnf[ic, ig]

        Jv += JvCnf[ic, ig]
        C[ic, ig] = copy.deepcopy(g[ic, ig].property('C'))
        if Jv >= 0.0: JCbase += JvCnf[ic, ig] * C[ic, ig][1]
    F += w_cnf * (Jv - Jexp[ic]) ** 2.0  # / S_g[ic]**2.0 #/ float(len(cut_n_flow_length) + 1) #/ Jexp[ic]**2.0
    if JCbase > 0.0: C_base = JCbase / Jv

    for ic in range(1, len(cut_n_flow_length) + 1):
        Jv = 0.0
        JCbase = 0.0
        data = row = col = None
        for ig in range(Nb_of_roots):
            g[ic, ig] = flux.flux(g[ic, ig], psi_e = psi_base + DP_cnf[ic], psi_base = psi_base, invert_model = True)
            g[ic, ig] = init_some_MTG_properties(g[ic, ig], tau = Js, Cini = Cini, Cpeg_ini = Cpeg_ini, t = 1, Ps = Ps)
            nb_v = g[ic, ig].nb_vertices()
            Fdx = 1.0
            Fdx_old = 1.
            Jv_old = 1.
            # Newton-Raphson schemes: in pressure_calculation calculation of dx,
            # array with dP, dC and dCpeg (if any) variation of the variables between two Newton step. Then the Newton scheme stops when
            # Fdx > eps see below
            while Fdx > eps:
                g[ic, ig], dx, data, row, col = routine_calculation(g[ic, ig], sigma = sigma,
                                                                    Ce = Ce, Pe = parameter.exp['psi_e'],
                                                                    Pbase = parameter.exp['psi_base'],
                                                                    Cse = Cse, dP = DP_cnf[ic])
                Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
                JvCnf[ic, ig] = g[ic, ig].property('J_out')[1]
                # if Flag_verbose: print local_j, Fdx, (Fdx - Fdx_old)
                if abs(JvCnf[ic, ig] - Jv_old) < 1.0e-4:
                    break
                if abs(Fdx - Fdx_old) < eps:
                    break
                Fdx_old = Fdx
                Jv_old = JvCnf[ic, ig]

            Jv += JvCnf[ic, ig]
            C[ic, ig] = copy.deepcopy(g[ic, ig].property('C'))
            if Jv >= 0.0: JCbase += JvCnf[ic, ig] * C[ic, ig][1]
        F += w_cnf * (Jv - Jexp[ic]) ** 2.0  # / S_g[ic]**2.0 #/ float(len(cut_n_flow_length) + 1) #/ Jexp[ic]**2.0
        if JCbase > 0.0: C_base = JCbase / Jv

    return JvCnf, F, C

###############################################################################
# Main
###############################################################################

if __name__ == '__main__':

    # read input parameters from results csv files

    # dir_csv = '/home/fabrice/Documents/hydroroot_FB-3/example/'
    # # parameter.solute['Sigma'] = 1.0
    # filename = filename.replace('parameters_', '')
    # fn = dir_csv + filename.replace('.yml', '')  + '.csv'
    # df_csv = pd.read_csv(fn, sep = ',', keep_default_na = True)
    # parameter.solute['J_s'] = df_csv['J_s'][0]
    # parameter.solute['P_s'] = df_csv['P_s'][0]
    # parameter.hydro['k0'] = df_csv['k'][0]

    # files with the experimental data of cnf and Jv(P) see DP_cnf, Jexp, list_DP and list_Jext
    # see comments at the top of this file
    fn = 'data/maize_cnf_data.csv'
    df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)
    fn = 'data/maize_Lpr_data.csv'
    df_exp2 = pd.read_csv(fn, sep = ',', keep_default_na = True)

    # files names of the architecture if reconstructed from a file
    # if not we just give a dummy name for the loop used to launch run
    filename = []
    parameter.archi['seed'] = [1]
    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))
    
    # dataframe used to save and export results: cnf and Jv(P)
    _col_names = ['max_length', 'Jexp cnf', 'Jv cnf', 'surface', 'length']
    results = {}
    for key in _col_names:
        results[key] = []
    _col_names2 = ['dp', 'Jexp(P)', 'Jv(P)', 'Cbase']
    results2 = {}
    for key in _col_names2:
        results2[key] = []

    ############################
    # get value from yaml file
    ############################

    # primary_length.append(parameter.archi['max_length'][0])
    delta = parameter.archi['branching_delay'][0]
    nude_length = parameter.archi['nude_length'][0]
    seed = parameter.archi['seed'][0]
    axfold = parameter.output['axfold'][0]
    radfold = parameter.output['radfold'][0]

    # Conductancies: mananging the fact there are or not different values between the primary and laterals
    # and the fact there are multiply by axfold and radfold
    k = []
    if type(parameter.hydro['k0']) != list:
        k.append(parameter.hydro['k0'] * radfold)
        k.append(parameter.hydro['k0'] * radfold)
    else:
        if len(parameter.hydro['k0'])>1:
            k.append(parameter.hydro['k0'][0] * radfold)
            k.append(parameter.hydro['k0'][1] * radfold)
        else:
            k.append(parameter.hydro['k0'][0] * radfold)
            k.append(parameter.hydro['k0'][0] * radfold)

    exp_axial = parameter.hydro['axial_conductance_data']
    axial_data = ([exp_axial[0], exp_axial[1]])
    axial_data = list(axial(axial_data, axfold))
    if len(exp_axial) == 4:
        axial_lr = ([exp_axial[2], exp_axial[3]])
        axial_lr = list(axial(axial_lr, axfold))
    else:
        axial_lr = copy.deepcopy(axial_data)

    J_s = parameter.solute['J_s']
    P_s = parameter.solute['P_s']
    Cse = parameter.solute['Cse'] * 1e-9 # mol/m3 -> mol/microL, external permeating solute concentration
    Ce = parameter.solute['Ce'] * 1e-9 # mol/m3 -> mol/microL, external non-permeating solute concentration
    Cini = Cse # initialization solute concentration into the xylem vessels
    Cpeg_ini = Ce # initialization non-permeating solute concentration into the xylem vessels: not 0.0 because more num instability
    Sigma = parameter.solute['Sigma'] # reflection coefficient, fixed in this script
    Pi_e_peg = osmotic_p_peg(Ce, unit_factor = 8.0e6)  # from Ce mol/microL to g/g, external osmotic pressure of non-permeating in MPa
    C_base_ini = 0.0

    data = None
    row = None
    col = None
    w_cnf = w_Lpr = 1. #weight on cnf cost function

    # functions that resolve the matrix system used in the Newton-Raphson scheme
    # different function depending on the presence of non-permeating solute, because there is one unknown less Cpeg
    routine_calculation = None
    if Ce <= 0.:
        # no non-permeating solute present
        routine_calculation = pressure_calculation_no_non_permeating_solutes
    else:
        routine_calculation = pressure_calculation

    # the objective function calculation to call depending on the data we fit
    if Flag_data_to_use == "cnf":
        fun_objective = fun_cnf_only
    elif Flag_data_to_use == "JvP":
        fun_objective = fun_JvP_only
    else:
        fun_objective = fun

    dK_constraint = 0.0 #-0.04 # dK/dx >= dK_constraint
    dK_constraint_max = 6.0e-2 # deprecated
    _tol = 5.0e-7 # does not have significant impact !!?? used in some minimize.optimize solver
    eps = 1.0e-9 # global: stop criterion for the Newton-Raphson loop in Jv_P_calculation and Jv_cnf_calculation

    # Parameter bounds
    Kbnds = (1.0e-10, np.inf) # axial conductance
    kbnds = (0.0, np.inf) # radial conductivity
    Jbnds = (1e-15, np.inf) # Js
    Pbnds = (1e-15, np.inf) # Ps

    psi_base = parameter.exp['psi_base']
    # default value for the pressure difference between the external medium and the base
    DP_cnf = []
    DP_cnf.append(parameter.exp['psi_e'] - psi_base)

    # variables used for the results output see end of script
    K = {}
    K['x pr'] = axial_data[0]
    K['K1st pr'] = axial_data[1]
    dK = pd.DataFrame(K, columns = ['x pr', 'K1st pr'])
    K = {}
    K['x lr'] = axial_lr[0]
    K['K1st lr'] = axial_lr[1]
    dK1st = pd.DataFrame(K, columns = ['x lr', 'K1st lr'])
    dK1st = pd.concat([dK,dK1st], axis = 1).fillna("")

    # we can give a list of architecture file names
    for f in filename:
        # architecture file to dataframe
        df = read_archi_data(f) if parameter.archi['read_architecture'] else None
        print(f.replace(glob.glob(parameter.archi['input_dir'])[0],"") , seed)
        index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")

        # read the data measurements from data base, cut-n-flow: flux, cut length and pressure difference
        for key in df_exp['arch']:
            if str(key).lower() in index.lower():
                _list = df_exp[df_exp.arch == key].filter(regex = '^J').dropna(axis = 1).values.tolist()
                Jexp = _list[0] # basal output flux
                _list = df_exp[df_exp.arch == key].filter(regex = '^lcut').dropna(axis = 1).values.tolist()
                cut_n_flow_length = _list[0] # cut lengthes
                _list = df_exp[df_exp.arch == key].filter(regex = '^dP').dropna(axis = 1).values.tolist()
                # the pressure difference is usually constant but sometimes, due to flow meter saturation, it may change
                # in that case a list of values is given
                if len(_list[0]) != 0:
                    DP_cnf = _list[0]
                    if len(DP_cnf) < len(cut_n_flow_length)+1: # if constant we create the list with the constant value
                        for i in range(1, len(cut_n_flow_length) + 1): DP_cnf.append(_list[0][0])

                parameter.exp['psi_e'] = psi_base + DP_cnf[0]

        # read the data measurements from data base Jv(P): flux, pressure
        for key in df_exp2['arch']:
            if str(key).lower() in index.lower():
                _list = df_exp2[df_exp2.arch == key].filter(regex = '^J').dropna(axis = 1).values.tolist()
                list_Jext = _list[0] # basal output flux
                _list = df_exp2[df_exp2.arch == key].filter(regex = '^dP').dropna(axis = 1).values.tolist()
                list_DP = _list[0] # delta pressure
                dlpr = pd.DataFrame(list(zip(list_DP, list_Jext)), columns = ['dP', 'Jv'])
                # dlpr = dlpr.sort_values('dP')[dlpr['dP']>0.05]
                list_DP = list(dlpr['dP'])
                list_Jext = list(dlpr['Jv'])

        if Flag_w_Lpr: w_Lpr = 1.0 / len(list_DP)
        if Flag_w_cnf: w_cnf = len(cut_n_flow_length)

        # building the MTG
        ###################
        Nb_of_roots = 2 if "-L" in index else 1 # sometimes thera are 2 roots for a given measurement with seminals
        primary_length = 0.
        _length = 0
        _surface = 0
        for ig in range(Nb_of_roots):
            if ig == 1:
                f2 = f.replace("-L", "-R")
                df = read_archi_data(f2) if parameter.archi['read_architecture'] else None

            g_cut[0, ig], _p, _l, _s, _seed = root_builder(primary_length, df = df, segment_length = parameter.archi['segment_length'],
                length_data = parameter.archi['length_data'],  branching_variability = parameter.archi['branching_variability'],
                order_max = parameter.archi['order_max'], order_decrease_factor = parameter.archi['order_decrease_factor'],
                ref_radius = parameter.archi['ref_radius'], Flag_radius = Flag_radius)
            if _p > primary_length: primary_length = _p
            _length += _l
            _surface += _s
            base = {}
            for v in g_cut[0, ig]:
                base[v] = next(axis(g_cut[0, ig], v))
            g_cut[0, ig].properties()['axisbase'] = base
            S_g.append(_s)

            # case where the primary is shorter than laterals
            max_length = primary_length
            mylength = g_cut[0, ig].property('mylength')
            if max(mylength.values()) > max_length: max_length = max(mylength.values())

            # set conductance
            g_cut[0, ig] = set_conductances(g_cut[0, ig], axial_pr = axial_data, k0_pr = k[0], axial_lr = axial_lr, k0_lr = k[1])
            # flux calculation without solute transport a way to initialize
            g_cut[0, ig] = flux.flux(g_cut[0, ig], psi_e = psi_base + DP_cnf[0], psi_base = psi_base, invert_model = True)

            # add properties specific to solute transport
            g_cut[0, ig].add_property('C') # permeating solute concentration
            g_cut[0, ig].add_property('Cpeg') # non-permeating solute concentration needed if cut-n-flow with them in the medium
            g_cut[0, ig].add_property('theta') # see init_some_MTG_properties
            g_cut[0, ig].add_property('J_s') # see init_some_MTG_properties, at a certain time I tried varying Js with C
            g_cut[0, ig].add_property('P_s') # see init_some_MTG_properties, at a certain time I tried varying Js with C
            g_cut[0, ig].add_property('original_vid') # the indices change between the full root and the cut root a way
                                                      # to retrieve the original index see set_K_and_k
            g_cut[0, ig].add_property('mu') # the viscosity of the sap because could change from the water value when
                                            # non-permeating solute enter the cut roots

            # a simple record of the original vertex number in the full architecture
            # do this because below when we cut we reindex because equations system is resolved in matrix form on the
            #  so the vertices need to have proper indices
            # MTG
            d = {k: k for k in g_cut[0, ig].vertices(g_cut[0, ig].max_scale())}
            g_cut[0, ig].properties()['original_vid'] = d
            # ############ longitudinal CUTS ####################################
            ic = 1
            for cut_length in cut_n_flow_length:
                # print(cut_length)
                tip_id[ic, ig] = \
                    flux.segments_at_length(g_cut[0, ig], cut_length, dl = parameter.archi['segment_length'])
                g_cut[ic, ig] = flux.cut(g_cut[0, ig], cut_length, parameter.archi['segment_length'])
                for i in tip_id[ic, ig]:
                    v = g_cut[0, ig].parent(i)
                    g_cut[ic, ig].property('label')[v] = 'cut' # labelling the vertices at cut ends

                # Below reindex because the system is resolved in matrix form on the MTG so the vertices need to have proper indices
                g_cut[ic, ig].reindex()
                i = 0
                tip_id[ic, ig] = [] # reinitializing because the cut can be at ramification then one parent for 2 different cut vertices
                for vid in g_cut[ic, ig].vertices_iter(g_cut[ic, ig].max_scale()):
                    if g_cut[ic, ig].label(vid) == 'cut':
                        tip_id[ic, ig].append(vid)
                        i += 1
                g_cut[ic, ig], surface = radius.compute_surface(g_cut[ic, ig])
                S_g.append(surface)
                ic += 1

        # Optimization
        ##############
        # the parameter are normalized with their inital values to limit scale effect between them, not the best the best would
        # be to write the equation in dimensionless form but historicaly hydroroot was not written this way
        iKpr = iKlr = iJs = iPs =ikpr = iklr = 0 # indices used to select the correct parameters in the array x, see fun for instance
        if Flag_Optim:
            # setting bounds and initial values
            ix = -1
            bnds = [] # list of tuple for bounds
            xini_list = [] # list of initial values of parameters
            if Flag_Optim_K:
                for var in axial_data[1]:
                    xini_list.append(var)
                ix += len(axial_data[1])  # be careful for axial_data and axial_lr the indices will be use as end of list interval selection => +1
            iKpr = int(ix)

            if Flag_Optim_Klr:
                for var in axial_lr[1]:
                    xini_list.append(var)
                ix += len(axial_lr[1])
                iKlr = int(ix)

            if Flag_Optim_K or Flag_Optim_Klr:
                for i, val in enumerate(xini_list):
                    bnds.append(Kbnds)

            if Flag_Optim_Js:
                xini_list.append( J_s)
                bnds.append(Jbnds)
                ix += 1
                iJs = int(ix)
            if Flag_Optim_Ps:
                xini_list.append(P_s)
                bnds.append(Pbnds)
                ix += 1
                iPs = int(ix)
            if Flag_Optim_k:
                xini_list.append(k[0])
                bnds.append(kbnds)
                ix += 1
                ikpr = int(ix)
            if Flag_Optim_klr:
                xini_list.append(k[1])
                bnds.append(kbnds)
                ix += 1
                iklr = int(ix)
            if Flag_Optim_sigma:
                xini_list.append(Sigma)
                if Sigma>0.0:
                    b = 1.0/Sigma
                else:
                    b = 1.0
                bnds.append((0.0, b))
                ix += 1
                isig = int(ix)

            xini = np.array(xini_list)
            x = np.ones(len(xini)) # the array of parameter that will be optimized, equal unity because we optimize the
                                   # the parameters normalized by their initial value

            # array used for constraints see optimize.minimize doc
            n = len(x)
            n1 = len(axial_data[1])
            # linear constraints lb <= A.dot(x) <= ub
            A = np.zeros((n, n))
            lb = np.full(n, -np.inf)
            ub = np.full(n, np.inf)
            l = copy.deepcopy(parameter.hydro['axial_conductance_data'][0])
            if Flag_Optim_Klr:
                l.append(0)
                l.append(0)
            if Flag_Optim_k: l.append(0)
            if Flag_Optim_Klr: l.append(0)

            if Flag_Optim_K and Flag_Constraint:
                a = dK_constraint # constraint on the 1st derivative
                for i in range(n1 - 1): # downward derivative
                    A[i, i] = -1.
                    A[i, i + 1] = 1.
                    lb[i] = a * (l[i+1]-l[i])
                    # ub[i] = dK_constraint_max * (l[i + 1] - l[i])
                ineq_cons = ({'type': 'ineq', 'fun': fun_constraint}) # !! works for K, k, Ps and Js optimized
            else:
                # for the COLBYLA solver bounds are not managed as other see fun_bound_cobyla
                ineq_cons = {'type': 'ineq', 'fun': fun_bound_cobyla}

            constraints = optimize.LinearConstraint(A, lb, ub) if Flag_Constraint else None

            # res = optimize.minimize(fun_objective, x, bounds = bnds, method = 'trust-constr', options={'finite_diff_rel_step': 1e-1})
            # res = optimize.minimize(fun_objective, x, bounds = bnds, method = 'SLSQP', constraints = [ineq_cons], options={'ftol': 1.0e-9, 'eps': 1e-1})
            res = optimize.minimize(fun_objective, x, method='COBYLA', constraints = [ineq_cons])
            # res = optimize.minimize(fun_objective, x, method='TNC', bounds = bnds)
            # res = optimize.minimize(fun_objective, x, bounds = bnds, options={'ftol': _tol, 'eps': 1e-1})
            # res = optimize.minimize(fun_objective, x, bounds = bnds, method='nelder-mead', options={'fatol': 1.0e-9})
            # res = optimize.basinhopping(fun_objective, x, minimizer_kwargs = {"bounds" : bnds, "constraints" : constraints, "options" : {'ftol': _tol, 'eps': 1e-1}})
            # res = optimize.dual_annealing(fun_objective, bounds = bnds)

            # optimization results to parameters
            n = len(axial_data[1])
            _x = res.x * xini
            if Flag_Optim_K:
                axial_data[1] = list(_x[:iKpr + 1])
            if not Flag_Optim_Klr:
                axial_lr = axial_data
            else:
                axial_lr[1] = _x[iKpr + 1:iKlr + 1]

            if Flag_Optim_Js: J_s = _x[iJs]
            if Flag_Optim_Ps: P_s = _x[iPs]

            if Flag_Optim_klr:
                k[1] = _x[iklr]
            else:
                k[1] = None
            if Flag_Optim_k:
                k[0] = _x[ikpr]

            if Flag_Optim_sigma:
                Sigma = _x[isig]

            print(res.x)

    # Direct simulation with the optimized values or the values from the yaml file if no optimization asked
    g_cut = set_K_and_k(g_cut, axial_data, k[0], axial_lr = axial_lr, k_lr = k[1], nr = Nb_of_roots,
                        nl = len(cut_n_flow_length))

    if Flag_data_to_use in ['all', 'JvP']:
        print('****** JvP ******')
        JvP, F_JvP, C = Jv_P_calculation(g_cut, Sigma, J_s, P_s)
        for idP in range(len(list_DP)):
            C_min = np.inf
            C_max = -np.inf
            Jv = 0.0
            JCbase = 0.0
            C_base = 0.0
            for ig in range(Nb_of_roots):
                if min(C[idP, ig].values()) < C_min: C_min = min(C[idP, ig].values())
                if max(C[idP, ig].values()) > C_max: C_max = max(C[idP, ig].values())
                if JvP[idP, ig] >= 0.0: JCbase += JvP[idP, ig] * C[idP, ig][1] # the solute flux at the outlet is the sum as the water flux
                Jv += JvP[idP, ig]
            if JCbase > 0.0: C_base = JCbase / Jv # the outlet solute concentration : (j0 * C0 + j1 * C1) / (j0 + j1)

            results2['dp'].append(list_DP[idP])
            results2['Jv(P)'].append(Jv)
            results2['Jexp(P)'].append(list_Jext[idP])
            results2['Cbase'].append(C_base * 1e9)
            print(list_DP[idP], JvP[idP, ig], 'Cmin: ', C_min * 1e9, 'Cmax: ', C_max * 1e9, 'Cbase: ', C_base * 1e9)

    if Flag_data_to_use in ['all', 'cnf']:
        print('****** cut-n-flow ******')
        # # in cut and flow the outlet flow rate is positive then C_base is no longer an input boundary
        JvCnf, F_cnf, C = Jv_cnf_calculation(g_cut, Sigma, J_s, P_s)
        for ic in range(len(cut_n_flow_length) + 1):
            _surface = 0.
            _length = 0.
            C_min = np.inf
            C_max = -np.inf
            Jv = 0.0
            JCbase = 0.0
            C_base = 0.0
            for ig in range(Nb_of_roots):
                g_cut[ic, ig], _s = radius.compute_surface(g_cut[ic, ig])
                _l = g_cut[ic, ig].nb_vertices(scale = 1) * parameter.archi['segment_length']
                _surface += _s
                _length += _l
                Jv += JvCnf[ic, ig]
                if min(C[ic, ig].values()) < C_min: C_min = min(C[ic, ig].values())
                if max(C[ic, ig].values()) > C_max: C_max = max(C[ic, ig].values())
                if JvCnf[ic, ig] >= 0.0: JCbase += JvCnf[ic, ig] * C[ic, ig][1] # the solute flux at the outlet is the sum as the water flux
            if JCbase > 0.0: C_base = JCbase / Jv

            if ic > 0: max_length = cut_n_flow_length[ic - 1]
            print(max_length, JvCnf[ic, ig], 'Cmin: ', C_min * 1e9, 'Cmax: ', C_max * 1e9, 'Cbase: ', C_base * 1e9)
            # results['plant'].append(index)
            results['max_length'].append(max_length)
            results['length'].append(_length)
            results['surface'].append(_surface)
            results['Jv cnf'].append(Jv)
            results['Jexp cnf'].append(Jexp[ic])

    ## just some parameter calculations for display
    # js_tot = 0
    # for ig in range(Nb_of_roots):
    #     g_cut[0, ig].add_property('DP')
    #     g_cut[0, ig].add_property('DC')
    #     g_cut[0, ig].add_property('jsurf')
    #     DP = g_cut[0, ig].property('DP')
    #     DC = g_cut[0, ig].property('DC')
    #     jsurf = g_cut[0, ig].property('jsurf')
    #     psi_in = g_cut[0, ig].property('psi_in')
    #     js = g_cut[0, ig].property('js')
    #     C = g_cut[0, ig].property('C')
    #     length = g_cut[0, ig].property('length')
    #     _radius = g_cut[0, ig].property('radius')
    #     j = g_cut[0, ig].property('j')
    #     for v in g_cut[0, ig].vertices_iter(scale = 1):
    #         js[v] = _radius[v] * 2 * np.pi * length[v] * (J_s + P_s * (Cse-C[v]) * 1e9)
    #         js_tot += js[v]
    #         DC[v] = -(Cse - C[v])*1e9
    #         DP[v] = parameter.exp['psi_base'] + DP_cnf[0] - psi_in[v]
    #         # psi_in[v] -= psi_base # just to put in relative pressure
    #         jsurf[v] = j[v] / (_radius[v] * 2 * np.pi * length[v])

    # some plots
    ############
    plt.ion()
    fig, axs = plt.subplots(2, 2)

    dr=pd.DataFrame()
    dr2=pd.DataFrame()
    F = F2 = 0.0
    if Flag_data_to_use in ['all', 'JvP']:
        dr2 = pd.DataFrame(results2, columns = _col_names2)
        dr2.sort_values(['dp'], inplace=True)
        dr2.plot.scatter(x = 'dp', y = 'Jexp(P)', ax = axs[0, 0])
        dr2.plot.line(x = 'dp', y = 'Jv(P)', ax = axs[0, 0])
        j = np.array(dr2.loc[:, ['Jv(P)', 'Jexp(P)']])
        F2 = w_Lpr * np.sum(np.diff(j)**2.0)
        axs[0, 0].set_ylim(j.min(),j.max())

    if Flag_data_to_use in ['all', 'cnf']:
        dr = pd.DataFrame(results, columns = _col_names)
        dr.plot.scatter(x = 'max_length', y = 'Jexp cnf', ax = axs[0, 1])
        dr.plot.line(x = 'max_length', y = 'Jv cnf', ax = axs[0, 1])
        j = np.array(dr.loc[:, ['Jv cnf', 'Jexp cnf']])
        F = w_cnf * np.sum(np.diff(j)**2.0)
        axs[0, 1].set_ylim(j.min(),j.max())

    d = pd.concat([dr,dr2], axis = 1).fillna("")

    X = {}
    X['kpr'] = [k[0]]
    if k[1] is None: k[1] = k[0]
    X['klr'] = [k[1]]
    X['Js'] = [J_s]
    X['Ps'] = [P_s]
    X['F cnf'] = [F]
    X['F Lpr'] = [F2]
    dX = pd.DataFrame(X, columns = ['kpr', 'klr', 'Js', 'Ps', 'F cnf','F Lpr'])
    d = pd.concat([d,dX], axis = 1).fillna("")

    K = {}
    K['x pr'] = axial_data[0]
    K['K pr'] = axial_data[1]
    dK = pd.DataFrame(K, columns = ['x pr', 'K pr'])
    d = pd.concat([d,dK], axis = 1).fillna("")
    dK1st.plot.scatter(x = 'x pr', y = 'K1st pr', ax = axs[1, 0])
    dK.plot.line(x = 'x pr', y = 'K pr', ax = axs[1, 0])
    axs[1, 0].set_ylim(min(dK1st['K1st pr'].min(),dK['K pr'].min()),max(dK1st['K1st pr'].max(),dK['K pr'].max()))

    fig.suptitle(index + 'WIP', fontsize = 16)

    K = {}
    K['x lr'] = axial_lr[0]
    K['K lr'] = axial_lr[1]
    dK = pd.DataFrame(K, columns = ['x lr', 'K lr'])
    d = pd.concat([d,dK], axis = 1).fillna("")
    d = pd.concat([d,dK1st], axis = 1).fillna("") # 1st guesses see at the beginning
    d.to_csv(output, index = False)

    fig.patch.set_facecolor('lightgrey')
    fig.tight_layout()

    print('F cnf, Lpr, total: ', F, F2, F+F2)
    print(index, ', ', '{:0.2f}'.format(k[0]), ', ', '{:0.2e}'.format(J_s), ', ', '{:0.2e}'.format(P_s), ', ', ', '.join('{:0.2e}'.format(i) for i in axial_data[1]))
    print('running time is ', time.time() - start_time)
    if Flag_radius: print('******* used radii from archi file ********')
