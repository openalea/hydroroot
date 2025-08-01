import glob
import copy
import math

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from openalea.mtg.algo import axis
from scipy import optimize, constants
from pathlib import Path

from openalea.hydroroot.display import plot
from openalea.hydroroot.read_file import read_archi_data
from openalea.hydroroot.main import root_builder, hydroroot_flow
from openalea.hydroroot import radius, flux, conductance
from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.conductance import set_conductances, axial
from openalea.hydroroot.water_solute_transport import pressure_calculation, pressure_calculation_no_non_permeating_solutes, \
    init_some_MTG_properties, osmotic_p_peg

def water_solute_model(parameter, df_archi =None, df_law =None,
                       df_cnf = None, df_JvP = None, Data_to_Optim = None, Flag_verbose = False,
                        data_to_use = 'all', output = None, optim_method = 'COBYLA', Flag_debug = False,
                        Flag_radius = True, Flag_Constraint = True, dK_constraint = -3.0e-2, Flag_w_Lpr = False,
                        Flag_w_cnf = False):
    """
    Perform direct simulations or parameters adjustment to fit data of Jv(P) and/or cut and flow experiments.
    Water and solute transport. **Works with constant radial conductivity.**

    :param parameter: Parameter - (see :func: Parameters)
    :param df_archi: DataFrame (None) - DataFrame with the architecture data (see below structure description)
    :param df_law: DataFrame list (None) - DataFrame with the length law data (see below structure description)
    :param df_cnf: DataFrame (None) - cut and flow data to fit (see below structure description)
    :param df_JvP: DataFrame (None) - Jv(P) data to fit (see below structure description)
    :param Data_to_Optim: string list (None) - list of parameters to adjust, if None perform direct simulation, empty list is equivalent to ['K', 'k', 'Js', 'Ps']
    :param Flag_verbose: boolean (False) - if True print intermediary results, optimization details, final simulation outputs, etc.
    :param data_to_use: string ('all') - data to fit either 'JvP' (Jv(P)), 'cnf' (cut and flow), or 'all' both
    :param output: string (None) - if not None output filename
    :param optim_method: string ('COBYLA') - solver method used in scipy.optimize.minimize see docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    :param Flag_debug: boolean (False) - if True print intermediary values of parameters adjustment
    :param Flag_radius: boolean (True) - if True use diameter recorded in architecture file if present, otherwise use ref_radius
    :param Flag_Constraint: boolean (True) - if True apply constraint on axial conductance 1st derivative (with COBYLA constraint may not be respected)
    :param dK_constraint: float (-3.0e-2) - lower bound of the axial conductance 1st derivative if Flag_Constraint = True (with COBYLA constraint may not be respected)
    :param Flag_w_Lpr: boolean (False)  - set weight to  1.0 / len(list_DP) on the JvP objective function
    :param Flag_w_cnf: boolean (False)  - set weight to  len(cut_n_flow_length) on the cut and flow objective function
    :return:
        - d: DataFrame with results
        - g_cut: dictionary with MTG at each cut

    **Data_to_Optim list**
        - 'K': optimize axial conductance K
        - 'k': optimize radial conductivity k
        - 'Js': optimize pumping rate Js
        - 'Ps': optimize permeability
        - 'Klr': optimize axial conductance of laterals Klr if <> than PR
        - 'klr': optimize radial conductivity of laterals klr if <> than PR
        - 'sigma': optimize the reflection coefficient

        ['K', 'k', 'Js', 'Ps'] will optimize the four parameters in the list.
        But, having too many parameters with this routine using local optimizer from scipy may not lead to any results.

    **df_archi column names**
        - distance_from_base_(mm), lateral_root_length_(mm), order

    **df_law**
        - list of 2 dataframe with the length law data: the first for the 1st order laterals on the primary root, the
          2nd for the laterals on laterals whatever their order (2nd, 3rd, ...)
        - column names: LR_length_mm , relative_distance_to_tip

    **df_cnf column names**
        - arch: sample name that must be contained in the 'input_file' of the yaml file
        - dP_Mpa: column with the working cut and flow pressure (in relative to the base) if constant, may be empty see below
        - J0, J1, ..., Jn: columns that start with 'J' containing the flux values, 1st the for the full root, then 1st cut, 2d cut, etc.
        - lcut1, ...., lcutn: columns starting with 'lcut' containing the maximum length to the base after each cut, 1st cut, 2d cut, etc. (not the for full root)
        - dP0, dP1,.., dPn: column starting with 'dP' containing the working pressure (in relative to the base) of each steps (if not constant): full root, 1st cut, 2d cut, etc.

    **df_JvP column names**
        - arch: sample name that must be contained in the 'input_file' of the yaml file
        - J0, J1, ..., Jn: columns that start with 'J' containing the flux values of each pressure steps
        - dP0, dP1,.., dPn: column starting with 'dP' containing the working pressure (in relative to the base) of each steps

    **outputfile (csv)**
        - column names: 'max_length', 'Jexp cnf (uL/s)', 'Jv cnf (uL/s)', 'surface (m2)', 'length (m)', 'dp', 'Jexp(P)', 'Jv(P)', 'Cbase', 'kpr', 'klr', 'Js', 'Ps', 'F cnf','F Lpr', 'x pr', 'K1st pr', 'x lr', 'K1st lr', 'K pr', 'K lr'

        i.e.: max length from the cut to the base, J cnf exp, J cnf sim, root surface, total root length, pressure (Jv(P),
        J exp Jv(P), J sim Jv(P), solute concentration at the base (Jv(P)), radial conductivity PR, radial conductivity LR,
        pumping rate, permeability, objective fct cnf, objective fct Jv(P), x pr distance to tip for PR, K 1st guess for PR,
        the same for laterals if different, then axial conductances for PR and LR.

    **Remark**
        - radial conductivity: single value or list of 2 values [kpr, klr]: 1st value for PR and 2nd for LR
        - axial conductance: it is possible to add K for LR
    """
    if Data_to_Optim is None:
        Data_to_Optim = []  # direct simulation
    elif len(Data_to_Optim) == 0:
        Data_to_Optim = ['K', 'k', 'Js', 'Ps']  # if Data_to_Optim = []

    Flag_Optim_K = ('K' in Data_to_Optim)  # optimize axial conductance K
    Flag_Optim_Klr = ('Klr' in Data_to_Optim)  # optimize axial conductance of laterals Klr if <> than PR
    Flag_Optim_klr = ('klr' in Data_to_Optim)  # optimize radial conductivity of laterals klr if <> than PR
    Flag_Optim_Js = ('Js' in Data_to_Optim)  # optimize pumping rate Js
    Flag_Optim_Ps = ('Ps' in Data_to_Optim)  # optimize permeability
    Flag_Optim_k = ('k' in Data_to_Optim)  # optimize radial conductivity k
    Flag_Optim_sigma = ('sigma' in Data_to_Optim)

    if df_cnf is None:
        # force on JvP data
        data_to_use = 'JvP'
    elif df_JvP is None:
        # force on cut and flow data
        data_to_use = 'cnf'
    # if both are none the simulation will be done using psi_base en ext given parameter and run as one JvP data point

    results = {}
    g_cut = {}
    tip_id = {}
    S_g = []
    cut_n_flow_length = []
    Jexp = []
    result_cv = []

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
            c[i] = x[i + 1] - x[i] - dK_constraint * (l[i + 1] - l[i])

        # bounds as inequality constraints needed by solver 'COBYLA' but redundant for 'SLSQP' with bounds
        for i in range(n, len(x)):
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
            axial_data[1] = list(
                _x[:iKpr + 1])  # np array : _x[:n] means the n 1st elements index [0;n-1] but _x[n] means _x index n

        if not Flag_Optim_Klr:
            _axial_lr = None
        else:
            _axial_lr = axial_lr
            _axial_lr[1] = _x[iKpr + 1:iKlr + 1]

        if Flag_Optim_Js:
            Js = _x[iJs]
        else:
            Js = J_s
        if Flag_Optim_Ps:
            Ps = _x[iPs]
        else:
            Ps = P_s

        if Flag_Optim_klr:
            klr = _x[iklr]
        else:
            klr = k[1]

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
        # n = len(axial_data[1])
        # c = np.ones(n - 1)
        # l = axial_data[0]
        # for i in range(n - 1):
        #     if _x[i+1] - _x[i] - dK_constraint * (l[i+1] - l[i]) < 0.0:
        #         F += 1.0

        if Flag_debug: print('{:0.3e}'.format(F), ' '.join('{:0.3e}'.format(i) for i in _x))

        result_cv.append([F, kpr, klr, Ps, Js] + axial_data[1])
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
            axial_data[1] = list(_x[:iKpr + 1])  # np array : _x[:n] means the n 1st elements index [0;n-1] but _x[n] means _x index n

        if not Flag_Optim_Klr:
            _axial_lr = None
        else:
            _axial_lr = axial_lr
            _axial_lr[1] = _x[iKpr + 1:iKlr + 1]

        if Flag_Optim_Js:
            Js = _x[iJs]
        else:
            Js = J_s

        if Flag_Optim_Ps:
            Ps = _x[iPs]
        else:
            Ps = P_s

        if Flag_Optim_klr:
            klr = _x[iklr]
        else:
            klr = k[1]

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
        # run simulation Jv(P)
        JvP, F, C = Jv_P_calculation(g, sigma, Js, Ps)

        if Flag_debug: print('{:0.3e}'.format(F), ' '.join('{:0.3e}'.format(i) for i in _x))

        result_cv.append([F, kpr, klr, Ps, Js] + axial_data[1])
        return F

    def fun_cnf_only(x):
        """
        Calculation of the objective function (Residual sum of squares) done on CnF data

        :param x: array of the parameters to optimize
        :return: F (float), the objective function
        """

        # Kx pr, Kx lr, Js, Ps,k lr and k pr may be optimized depends on the Flags_Optim_
        _x = x * xini
        if iKpr > 0: axial_data[1] = list(
            _x[:iKpr + 1])  # np array : _x[:n] means the n 1st elements index [0;n-1] but _x[n] means _x index n

        if not Flag_Optim_Klr:
            _axial_lr = None
        else:
            _axial_lr = axial_lr
            _axial_lr[1] = _x[iKpr + 1:iKlr + 1]

        if Flag_Optim_Js:
            Js = _x[iJs]
        else:
            Js = J_s
        if Flag_Optim_Ps:
            Ps = _x[iPs]
        else:
            Ps = P_s

        if Flag_Optim_klr:
            klr = _x[iklr]
        else:
            klr = k[1]

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
        # run simulation Jv(P)
        JvCnf, F, C = Jv_cnf_calculation(g, sigma, Js, Ps)

        if Flag_debug: print('{:0.3e}'.format(F), ' '.join('{:0.3e}'.format(i) for i in _x))

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
                    g[ic, ig].property('K_exp')[v] = g[0, ig].property('K_exp')[
                        vid]  # needed in pressure_calculation if Cpeg because calculation of K
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
        data = None
        row = None
        col = None
        for idP in range(len(list_DP)):
            Jv = 0.0
            for ig in range(Nb_of_roots):
                g[0, ig] = flux.flux(g[0, ig], psi_e = psi_base + list_DP[idP], psi_base = psi_base,
                                     invert_model = True)
                g[0, ig] = init_some_MTG_properties(g[0, ig], tau = Js, Cini = Cini, t = 1, Ps = Ps)
                nb_v = g[0, ig].nb_vertices()
                Fdx = 1.0
                Fdx_old = 1.
                Jv_old = 1.
                # Newton-Raphson schemes: in pressure_calculation_no_non_permeating_solutes calculation of dx,
                # array with dP and dC variation of the variables between two Newton step. Then the Newton scheme stops when
                # Fdx > eps see below
                while Fdx > eps:
                    g[0, ig], dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g[0, ig],
                                                                                                  sigma = sigma,
                                                                                                  Ce = Ce,
                                                                                                  Pe = parameter.exp[
                                                                                                      'psi_e'],
                                                                                                  Pbase = parameter.exp[
                                                                                                      'psi_base'],
                                                                                                  Cse = Cse,
                                                                                                  dP = list_DP[idP],
                                                                                                  C_base = None)
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
            F += w_Lpr * (Jv - list_Jext[idP]) ** 2.0

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
                                                                                               Pe = parameter.exp[
                                                                                                   'psi_e'],
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
        F += w_cnf * (Jv - Jexp[ic]) ** 2.0

        for ic in range(1, len(cut_n_flow_length) + 1):
            Jv = 0.0
            data = row = col = None
            for ig in range(Nb_of_roots):
                g[ic, ig] = flux.flux(g[ic, ig], psi_e = psi_base + DP_cnf[ic], psi_base = psi_base,
                                      invert_model = True)
                g[ic, ig] = init_some_MTG_properties(g[ic, ig], tau = Js, Cini = Cini, Cpeg_ini = Cpeg_ini, t = 1,
                                                     Ps = Ps)
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
                    # if Flag_debug: print local_j, Fdx, (Fdx - Fdx_old)
                    if abs(JvCnf[ic, ig] - Jv_old) < 1.0e-4:
                        break
                    if abs(Fdx - Fdx_old) < eps:
                        break
                    Fdx_old = Fdx
                    Jv_old = JvCnf[ic, ig]

                Jv += JvCnf[ic, ig]
                C[ic, ig] = copy.deepcopy(g[ic, ig].property('C'))
            F += w_cnf * (Jv - Jexp[ic]) ** 2.0

        return JvCnf, F, C


    # architecture file to dataframe
    index=''
    # if (df_archi is None) and parameter.archi['read_architecture']:
    #     # architecture with filename in aqua team format
    #     archi_f = glob.glob(parameter.archi['input_dir'] + parameter.archi['input_file'][0])
    #     archi_f = archi_f[0]
    #     df_archi = read_archi_data(archi_f) if parameter.archi['read_architecture'] else None
    #     index = archi_f.replace(glob.glob(parameter.archi['input_dir'])[0], "")
    
    if parameter.archi['read_architecture']:
        # architecture with filename in aqua team format
        # archi_f = glob.glob(parameter.archi['input_dir'] + parameter.archi['input_file'][0])
        # archi_f = archi_f[0]
        fname = str(Path(parameter.archi['input_dir']) / parameter.archi['input_file'][0]) # deals with path in windows
        archi_f = glob.glob(fname) # deals with wildcards ? deprecated ???
        archi_f = archi_f[0]
        
        if df_archi is None: df_archi = read_archi_data(archi_f)
        index = archi_f.replace(glob.glob(parameter.archi['input_dir'])[0], "")
    index = parameter.archi['input_file'][0]

    # length law data: override if necessary
    if df_law is not None:
        parameter.archi['length_data'] = df_law

    # dataframe used to save and export results: cnf and Jv(P)
    _col_names = ['max_length', 'Jexp cnf (uL/s)', 'Jv cnf (uL/s)', 'surface (m2)', 'length (m)']
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
        k.append(None)
    else:
        if len(parameter.hydro['k0']) > 1:
            k.append(parameter.hydro['k0'][0] * radfold)
            k.append(parameter.hydro['k0'][1] * radfold)
        else:
            k.append(parameter.hydro['k0'][0] * radfold)
            k.append(None)

    exp_axial = parameter.hydro['axial_conductance_data']
    axial_data = ([exp_axial[0], exp_axial[1]])
    axial_data = list(axial(axial_data, axfold))
    if len(exp_axial) == 4:
        axial_lr = ([exp_axial[2], exp_axial[3]])
        axial_lr = list(axial(axial_lr, axfold))
    else:
        axial_lr = None #copy.deepcopy(axial_data)

    J_s = parameter.solute['J_s']
    P_s = parameter.solute['P_s']
    Cse = parameter.solute['Cse'] * 1e-9  # mol/m3 -> mol/microL, external permeating solute concentration
    Ce = parameter.solute['Ce'] * 1e-9  # mol/m3 -> mol/microL, external non-permeating solute concentration
    Cini = Cse  # initialization solute concentration into the xylem vessels
    Cpeg_ini = Ce  # initialization non-permeating solute concentration into the xylem vessels: not 0.0 because more num instability
    Sigma = parameter.solute['Sigma']  # reflection coefficient, fixed in this script
    Pi_e_peg = osmotic_p_peg(Ce, unit_factor = 8.0e6)  # from Ce mol/microL to g/g, external osmotic pressure of non-permeating in MPa

    data = None
    row = None
    col = None
    w_cnf = w_Lpr = 1.  # weight on cnf cost function

    # functions that resolve the matrix system used in the Newton-Raphson scheme
    # different function depending on the presence of non-permeating solute, because there is one unknown less Cpeg
    routine_calculation = None
    if Ce <= 0.:
        # no non-permeating solute present
        routine_calculation = pressure_calculation_no_non_permeating_solutes
    else:
        routine_calculation = pressure_calculation

    # the objective function calculation to call depending on the data we fit
    if data_to_use == "cnf":
        fun_objective = fun_cnf_only
    elif data_to_use == "JvP":
        fun_objective = fun_JvP_only
    else:
        fun_objective = fun

    dK_constraint_max = 6.0e-2  # deprecated
    _tol = 5.0e-7  # does not have significant impact !!?? used in some minimize.optimize solver
    eps = 1.0e-9  # global: stop criterion for the Newton-Raphson loop in Jv_P_calculation and Jv_cnf_calculation

    # Parameter bounds
    Kbnds = (1.0e-10, np.inf)  # axial conductance
    kbnds = (0.01, np.inf)  # radial conductivity
    Jbnds = (1e-15, np.inf)  # Js
    Pbnds = (1e-15, np.inf)  # Ps

    psi_base = parameter.exp['psi_base']
    # default value for the pressure difference between the external medium and the base
    DP_cnf = []
    DP_cnf.append(parameter.exp['psi_e'] - psi_base)

    # variables used for the results output see end of script
    K = {}
    K['x pr'] = axial_data[0]
    K['K1st pr'] = axial_data[1]
    dK1st = pd.DataFrame(K, columns = ['x pr', 'K1st pr'])
    K = {}
    if axial_lr:
        K['x lr'] = axial_lr[0]
        K['K1st lr'] = axial_lr[1]
    else:
        K['x lr'] = axial_data[0]
        K['K1st lr'] = axial_data[1]
    dKlr1st = pd.DataFrame(K, columns = ['x lr', 'K1st lr'])

    # # architecture file to dataframe
    # df_archi = read_archi_data(archi_f) if parameter.archi['read_architecture'] else None
    # index = archi_f.replace(glob.glob(parameter.archi['input_dir'])[0], "")

    # read the data measurements from data base, cut-n-flow: flux, cut length and pressure difference
    if df_cnf is not None:
        for key in df_cnf['arch']:
            if str(key).lower() in index.lower():
                _list = df_cnf[df_cnf.arch == key].filter(regex = '^J').dropna(axis = 1).values.tolist()
                Jexp = _list[0]  # basal output flux
                _list = df_cnf[df_cnf.arch == key].filter(regex = '^lcut').dropna(axis = 1).values.tolist()
                cut_n_flow_length = _list[0]  # cut lengthes
                _list = df_cnf[df_cnf.arch == key].filter(regex = '^dP').dropna(axis = 1).values.tolist()
                # the pressure difference is usually constant but sometimes, due to flow meter saturation, it may change
                # in that case a list of values is given
                if len(_list[0]) != 0:
                    DP_cnf = _list[0]
                else:
                    DP_cnf = []
                    DP_cnf.append(parameter.exp['psi_e'] - psi_base) # for compatibility reason with first analysis on arabidopsis

                if len(DP_cnf) < len(cut_n_flow_length) + 1:  # if constant we create the list with the constant value
                    for i in range(1, len(cut_n_flow_length) + 1): DP_cnf.append(DP_cnf[0])

                parameter.exp['psi_e'] = psi_base + DP_cnf[0]

    # read the data measurements from data base Jv(P): flux, pressure
    if df_JvP is not None:
        for key in df_JvP['arch']:
            if str(key).lower() in index.lower():
                _list = df_JvP[df_JvP.arch == key].filter(regex = '^J').dropna(axis = 1).values.tolist()
                list_Jext = _list[0]  # basal output flux
                _list = df_JvP[df_JvP.arch == key].filter(regex = '^dP').dropna(axis = 1).values.tolist()
                list_DP = _list[0]  # delta pressure

                ## below juste to get data above a minimum dP
                # dlpr = pd.DataFrame(list(zip(list_DP, list_Jext)), columns = ['dP', 'Jv'])
                # dlpr = dlpr.sort_values('dP')[dlpr['dP']>0.05]
                # list_DP = list(dlpr['dP'])
                # list_Jext = list(dlpr['Jv'])
    else:
        list_Jext = [parameter.exp['Jv']]
        list_DP = [parameter.exp['psi_e'] - psi_base]

    if Flag_w_Lpr: w_Lpr = 1.0 / len(list_DP)
    if Flag_w_cnf: w_cnf = len(cut_n_flow_length)

    # building the MTG
    ###################
    Nb_of_roots = 2 if "-L" in index else 1  # sometimes thera are 2 roots for a given measurement with seminals
    if df_archi is None:
        primary_length = parameter.archi['primary_length'][0]
    else:
        primary_length = 0.
    _length = 0
    _surface = 0
    for ig in range(Nb_of_roots):
        if ig == 1:
            f2 = archi_f.replace("-L", "-R")
            df_archi = read_archi_data(f2) if parameter.archi['read_architecture'] else None

        g_cut[0, ig], _p, _l, _s, _seed = root_builder(df = df_archi,
                                                        primary_length = parameter.archi['primary_length'][0],
                                                        seed = parameter.archi['seed'][0],
                                                        delta = parameter.archi['branching_delay'][0],
                                                        nude_length = parameter.archi['nude_length'][0],
                                                        segment_length = parameter.archi['segment_length'],
                                                        length_data = parameter.archi['length_data'],
                                                        branching_variability = parameter.archi['branching_variability'],
                                                        order_max = parameter.archi['order_max'],
                                                        order_decrease_factor = parameter.archi['order_decrease_factor'],
                                                        ref_radius = parameter.archi['ref_radius'],
                                                        Flag_radius = Flag_radius)
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
        g_cut[0, ig] = set_conductances(g_cut[0, ig], axial_pr = axial_data, k0_pr = k[0], axial_lr = axial_lr,
                                        k0_lr = k[1])
        # flux calculation without solute transport a way to initialize
        g_cut[0, ig] = flux.flux(g_cut[0, ig], psi_e = psi_base + DP_cnf[0], psi_base = psi_base, invert_model = True)

        # add properties specific to solute transport
        g_cut[0, ig].add_property('C')  # permeating solute concentration
        g_cut[0, ig].add_property('Cpeg')  # non-permeating solute concentration needed if cut-n-flow with them in the medium
        g_cut[0, ig].add_property('theta')  # see init_some_MTG_properties
        g_cut[0, ig].add_property('J_s')  # see init_some_MTG_properties, at a certain time I tried varying Js with C
        g_cut[0, ig].add_property('P_s')  # see init_some_MTG_properties, at a certain time I tried varying Js with C
        g_cut[0, ig].add_property('original_vid')  # the indices change between the full root and the cut root a way
        # to retrieve the original index see set_K_and_k
        g_cut[0, ig].add_property('mu')  # the viscosity of the sap because could change from the water value when
        # non-permeating solute enter the cut roots

        # a simple record of the original vertex number in the full architecture
        # do this because below when we cut we reindex because equations system is resolved in matrix form on the
        #  so the vertices need to have proper indices
        # MTG
        d = {vid: vid for vid in g_cut[0, ig].vertices(g_cut[0, ig].max_scale())}
        g_cut[0, ig].properties()['original_vid'] = d
        # ############ longitudinal CUTS ####################################
        ic = 1
        for cut_length in cut_n_flow_length:
            # print(cut_length)
            tip_id[ic, ig] = flux.segments_at_length(g_cut[0, ig], cut_length, dl = parameter.archi['segment_length'])
            g_cut[ic, ig] = flux.cut(g_cut[0, ig], cut_length, parameter.archi['segment_length'])
            for i in tip_id[ic, ig]:
                v = g_cut[0, ig].parent(i)
                g_cut[ic, ig].property('label')[v] = 'cut'  # labelling the vertices at cut ends

            # Below reindex because the system is resolved in matrix form on the MTG so the vertices need to have proper indices
            g_cut[ic, ig].reindex()
            i = 0
            tip_id[ic, ig] = []  # reinitializing because the cut can be at ramification then one parent for 2 different cut vertices
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
    iKpr = iKlr = iJs = iPs = ikpr = iklr = 0  # indices used to select the correct parameters in the array x, see fun for instance
    if Data_to_Optim:
        # setting bounds and initial values
        ix = -1
        bnds = []  # list of tuple for bounds
        xini_list = []  # list of initial values of parameters
        if Flag_Optim_K:
            for var in axial_data[1]:
                xini_list.append(var)
            ix += len(axial_data[1])  # be careful with axial_data and axial_lr the indices will be use as end of list interval selection => +1
        iKpr = int(ix)

        if Flag_Optim_Klr:
            if axial_lr is None: axial_lr = copy.deepcopy(axial_data)
            for var in axial_lr[1]:
                xini_list.append(var)
            ix += len(axial_lr[1])
            iKlr = int(ix)

        if Flag_Optim_K or Flag_Optim_Klr:
            for i, val in enumerate(xini_list):
                bnds.append(Kbnds)

        if Flag_Optim_Js:
            xini_list.append(J_s)
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
            if k[1] is None: k[1] = copy.deepcopy(k[0])
            xini_list.append(k[1])
            bnds.append(kbnds)
            ix += 1
            iklr = int(ix)
        if Flag_Optim_sigma:
            xini_list.append(Sigma)
            if Sigma > 0.0:
                b = 1.0 / Sigma
            else:
                b = 1.0
            bnds.append((0.0, b))
            ix += 1
            isig = int(ix)

        xini = np.array(xini_list)
        x = np.ones(len(xini))  # the array of parameter that will be optimized, equal unity because we optimize the
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
            a = dK_constraint  # constraint on the 1st derivative
            for i in range(n1 - 1):  # downward derivative
                A[i, i] = -1.
                A[i, i + 1] = 1.
                lb[i] = a * (l[i + 1] - l[i])
                # ub[i] = dK_constraint_max * (l[i + 1] - l[i])
            ineq_cons = ({'type': 'ineq', 'fun': fun_constraint})  # !! works for K, k, Ps and Js optimized
        else:
            # for the COLBYLA solver bounds are not managed as other see fun_bound_cobyla
            ineq_cons = {'type': 'ineq', 'fun': fun_bound_cobyla}

        constraints = optimize.LinearConstraint(A, lb, ub) if Flag_Constraint else None

        if optim_method == 'COBYLA':
            res = optimize.minimize(fun_objective, x, method = optim_method, constraints = [ineq_cons])
        elif optim_method == 'SLSQP':
            res = optimize.minimize(fun_objective, x, bounds = bnds, method = optim_method, constraints = [ineq_cons],
                                    options = {'ftol': 1.0e-9, 'eps': 1e-1})
        else:
            res = optimize.minimize(fun_objective, x, bounds = bnds, method = optim_method)

        # res = optimize.minimize(fun_objective, x, bounds = bnds, method = 'trust-constr', options={'finite_diff_rel_step': 1e-1})
        # res = optimize.minimize(fun_objective, x, method='TNC', bounds = bnds)
        # res = optimize.minimize(fun_objective, x, bounds = bnds, options={'ftol': _tol, 'eps': 1e-1})
        # res = optimize.minimize(fun_objective, x, bounds = bnds, method='nelder-mead', options={'fatol': 1.0e-9})

        # optimization results to parameters
        n = len(axial_data[1])
        _x = res.x * xini
        if Flag_Optim_K:
            axial_data[1] = list(_x[:iKpr + 1])
        if Flag_Optim_Klr:
            axial_lr[1] = _x[iKpr + 1:iKlr + 1]

        if Flag_Optim_Js: J_s = _x[iJs]
        if Flag_Optim_Ps: P_s = _x[iPs]

        if Flag_Optim_klr:
            k[1] = _x[iklr]
            # else:
            #     k[1] = None
        if Flag_Optim_k:
            k[0] = _x[ikpr]

        if Flag_Optim_sigma:
            Sigma = _x[isig]

        if Flag_verbose: print(res.x)

    # Direct simulation with the optimized values or the values from the yaml file if no optimization asked
    g_cut = set_K_and_k(g_cut, axial_data, k[0], axial_lr = axial_lr, k_lr = k[1], nr = Nb_of_roots,
                        nl = len(cut_n_flow_length))

    if data_to_use in ['all', 'JvP']:
        JvP, F_JvP, C = Jv_P_calculation(g_cut, Sigma, J_s, P_s)
        for idP in range(len(list_DP)):
            Jv = 0.0
            C_base = 0.0
            for ig in range(Nb_of_roots):
                # C_base here is in the middle of the 1st MTG element because the boundary condition chosen here is
                # dC/dx = 0, so the concentration at the root boundary is the same.
                # if there are several roots (as when the experiment is done with 2 seminals), since the MTG elements
                # are equals the tital C_base is the average
                C_base += C[idP, ig][1] / float(Nb_of_roots)
                Jv += JvP[idP, ig]

            results2['dp'].append(list_DP[idP])
            results2['Jv(P)'].append(Jv)
            results2['Jexp(P)'].append(list_Jext[idP])
            results2['Cbase'].append(C_base * 1e9)

    if data_to_use in ['all', 'cnf']:
        JvCnf, F_cnf, C = Jv_cnf_calculation(g_cut, Sigma, J_s, P_s)
        for ic in range(len(cut_n_flow_length) + 1):
            _surface = 0.
            _length = 0.
            Jv = 0.0
            C_base = 0.0
            for ig in range(Nb_of_roots):
                g_cut[ic, ig], _s = radius.compute_surface(g_cut[ic, ig])
                _l = g_cut[ic, ig].nb_vertices(scale = 1) * parameter.archi['segment_length']
                _surface += _s
                _length += _l
                Jv += JvCnf[ic, ig]
                C_base += C[ic, ig][1] / float(Nb_of_roots) # C_base calculation see above JvP case

            if ic > 0: max_length = cut_n_flow_length[ic - 1]

            results['max_length'].append(max_length)
            results['length (m)'].append(_length)
            results['surface (m2)'].append(_surface)
            results['Jv cnf (uL/s)'].append(Jv)
            results['Jexp cnf (uL/s)'].append(Jexp[ic])

    ## just some parameter calculations for display
    for ic in range(int(len(g_cut)/Nb_of_roots)): # FB 250328: added loop over cut
        js_tot = 0
        for ig in range(Nb_of_roots):
            g_cut[ic, ig].add_property('DP')
            g_cut[ic, ig].add_property('DC')
            g_cut[ic, ig].add_property('jsurf')
            g_cut[ic, ig].add_property('js')
            g_cut[ic, ig].add_property('DPi')
            DP = g_cut[ic, ig].property('DP')
            DC = g_cut[ic, ig].property('DC')
            DPi = g_cut[ic, ig].property('DPi')
            jsurf = g_cut[ic, ig].property('jsurf')
            psi_in = g_cut[ic, ig].property('psi_in')
            js = g_cut[ic, ig].property('js')
            C = g_cut[ic, ig].property('C')
            length = g_cut[ic, ig].property('length')
            _radius = g_cut[ic, ig].property('radius')
            j = g_cut[ic, ig].property('j')
            for v in g_cut[ic, ig].vertices_iter(scale = 1):
                js[v] = _radius[v] * 2 * np.pi * length[v] * (J_s + P_s * (Cse-C[v]) * 1e9)
                js_tot += js[v]
                DC[v] = -(Cse - C[v])*1e9
                DPi[v] = DC[v] * constants.R * 293 * 1e-6 # 1e-6 because in MPa
                DP[v] = parameter.exp['psi_base'] + DP_cnf[0] - psi_in[v]
                # psi_in[v] -= psi_base # just to put in relative pressure
                jsurf[v] = j[v] / (_radius[v] * 2 * np.pi * length[v])

    dr = pd.DataFrame()
    dr2 = pd.DataFrame()
    F = F2 = 0.0

    if Flag_verbose:
        pd.set_option('display.max_columns', None)
        pd.set_option('display.expand_frame_repr', False)

    if data_to_use in ['all', 'JvP']:
        dr2 = pd.DataFrame(results2, columns = _col_names2)
        # dr2.sort_values(['dp'], inplace=True)
        j = np.array(dr2.loc[:, ['Jv(P)', 'Jexp(P)']])
        F2 = w_Lpr * np.sum(np.diff(j) ** 2.0)
        if Flag_verbose:
            print('****** JvP ******')
            print(dr2)

    if data_to_use in ['all', 'cnf']:
        dr = pd.DataFrame(results, columns = _col_names)
        j = np.array(dr.loc[:, ['Jv cnf (uL/s)', 'Jexp cnf (uL/s)']])
        F = w_cnf * np.sum(np.diff(j) ** 2.0)
        if Flag_verbose:
            print('****** cut-n-flow ******')
            print(dr)


    d = pd.concat([dr, dr2], axis = 1).fillna("")

    X = {}
    X['kpr'] = [k[0]]
    if k[1] is None:
        X['klr'] = [k[0]]
    else:
        X['klr'] = [k[1]]
    X['Js'] = [J_s]
    X['Ps'] = [P_s]
    X['F cnf'] = [F]
    X['F Lpr'] = [F2]
    dX = pd.DataFrame(X, columns = ['kpr', 'klr', 'Js', 'Ps', 'F cnf', 'F Lpr'])

    K = {}
    # K['x pr'] = axial_data[0]
    K['K pr'] = axial_data[1]
    dK = pd.DataFrame(K, columns = ['K pr'])

    Klr = {}
    if axial_lr:
        Klr['K lr'] = axial_lr[1]
    else:
        Klr['K lr'] = axial_data[1]
    dKlr = pd.DataFrame(Klr, columns = ['K lr'])

    d = pd.concat([d, dX, dK1st, dK, dKlr1st, dKlr], axis = 1).fillna("")

    if output is not None: d.to_csv(output, index = False)

    if Flag_verbose:
        print('****** End ******')
        print('objective functions: ', 'F cnf: {:0.2e}'.format(F), 'F JvP: {:0.2e}'.format(F2), 'F tot: {:0.2e}'.format(F+F2))
        print(index, ',', 'k: {:0.2f}'.format(k[0]), ',', 'Js: {:0.2e}'.format(J_s), ',', 'Ps: {:0.2e}'.format(P_s), ', K: [',
              ', '.join('{:0.2e}'.format(i) for i in axial_data[1]), ']')

    return d, g_cut

def pure_hydraulic_model(parameter = Parameters(), df_archi = None, df_law =None, df_exp = None,
                         Data_to_Optim = None, output = None, Flag_verbose = False,
                        Flag_radius = False, Flag_Constraint = True, dK_constraint = -3.0e-2, dk_max = 0.1):
    """
    Perform direct simulations or parameters adjustment to fit data of cut and flow experiment.
    Water transport only, electrical network analogy

    :param parameter: Parameter - (see :func: Parameters)
    :param df_archi: DataFrame (None) - DataFrame with the architecture data (see below structure description)
    :param df_law: DataFrame list (None) - DataFrame with the length law data (see below structure description)
    :param df_exp: DataFrame (None) - data to fit
    :param Data_to_Optim: string list (None) - list of parameters to adjust, if None perform direct simulation, ['K', 'k']
    :param output: string (None) - if not None output filename
    :param Flag_verbose: boolean (False) - if True print intermediary results
    :param Flag_radius: boolean (False) - if True use diameter recorded in architecture file if present, otherwise use ref_radius
    :param Flag_Constraint: boolean (True) - if True apply constraint on axial conductance 1st derivative
    :param dK_constraint: float (-3.0e-2) - lower bound of the axial conductance 1st derivative if Flag_Constraint = True
    :param dk_max: float (0.1) - the convergence criteria on
    :return:
        - df: DataFrame with results
        - g_cut: dictionary with MTG at each cut


    **df_archi column names**
        - distance_from_base_(mm), lateral_root_length_(mm), order

    **df_law**
        - list of 2 dataframe with the length law data: the first for the 1st order laterals on the primary root, the
          2nd for the laterals on laterals whatever their order (2nd, 3rd, ...)
        - column names: LR_length_mm , relative_distance_to_tip

    **The adjustment is performed as follows**
        1. pre-optimization with the adjustment of axfold and radfold, K and k factor, if only k adjustment is asked then step 1 is not performed
        2. loop of two successive adjustments: 1st K adjustment then k adjustment. The loop stop when change of k is below dk_max

    **Data_to_Optim list of string**
        - 'K': optimize axial conductance K only
        - 'k': optimize radial conductivity k only
        - [ ] empty list or ['K', 'k']: optimize K and k

    **df_exp column names**
        - arch: sample name that must be contained in the 'input_file' of the yaml file
        - J0, ..., Jn: columns containing the flux values of the full root, 1st cut, 2d cut, etc.
        - lcut1, ...., lcutn: columns containing the maximum length to the base after each cut, 1st cut, 2d cut, etc. (the primary length of the full root is calculated from the architecture)


    **outputfile**
        - column names: 'plant', 'cut length (m)', 'max_length', 'k (10-8 m/s/MPa)', 'length (m)', 'surface (m2)', 'Jv (uL/s)', 'Jexp (uL/s)'
        - if 'K' in Data_to_Optim add the following: 'x', 'K 1st', 'K optimized' the initial and adjusted K(x)

    **Remark**
        The routine is designed to work with a single value (float) for parameter.hydro['k0'].

    **example**

    >>> parameter = Parameters()
    >>> filename='parameters_fig-2-B.yml'
    >>> parameter.read_file(filename)
    >>> fn = 'data/arabido_cnf_data.csv'
    >>> df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)
    >>> df = pure_hydraulic_model(parameter,df_exp=df_exp, Flag_verbose=True, Data_to_Optim = ['k', 'K'])

    """
    if Data_to_Optim is None:
        Data_to_Optim = []  # direct simulation
    elif len(Data_to_Optim) == 0:
        Data_to_Optim = ['K', 'k']  # if Data_to_Optim = []

    Flag_Optim_K = ('K' in Data_to_Optim)  # optimize axial conductance K
    Flag_Optim_k = ('k' in Data_to_Optim)  # optimize radial conductivity k

    g_cut = {}
    tip_id = {}
    cut_n_flow_length = []
    _tol = 1.0e-9

    def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None, psi_base = 0.1, psi_e = 0.1):
        if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
        if k_radial is None: k_radial = parameter.hydro['k0']
        # compute axial & radial
        Kexp_axial_data = conductance.axial(axial_data, axfold)
        k_radial_data = conductance.radial(k_radial, axial_data, radfold)

        ## Step function
        # k_radial_data = conductance.radial_step(k_radial,3.0,x_step = 0.02, dx = parameter.archi['segment_length'], scale = radfold)

        # compute local jv and psi, global Jv, Keq
        g, Keq, Jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                           k0 = k_radial,
                                           Jv = _Jv[0],
                                           psi_e = psi_e,
                                           psi_base = psi_base,
                                           axial_conductivity_data = Kexp_axial_data,
                                           radial_conductivity_data = k_radial_data)

        return g, Keq, Jv

    def fun1(x):
        """
        Simulation of the flux at the different cut lengths according to the new parameter value

        Implementation 1: only axfold (Kx factor) and radfold (k radial factor) are changed

        :param x: the array of adjusted parameters
        :return: F the sum((Jv - Jv_exp) ** 2.0)
        """
        axfold = x[0]
        radfold = x[-1]

        g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], radfold = radfold, axfold = axfold, psi_base = psi_base,
                                                  psi_e = psi_base + DP_cnf[0])
        F = (Jv - _Jv[0]) ** 2.0
        count = 1
        for cut_length in cut_n_flow_length:
            # _g = g_cut[str(cut_length)].copy() # not necessary and time-consuming

            for vid in g_cut[str(cut_length)].vertices_iter(g_cut['tot'].max_scale()):
                g_cut[str(cut_length)].property('K')[vid] = g_cut['tot'].property('K')[vid]
                g_cut[str(cut_length)].property('k')[vid] = g_cut['tot'].property('k')[vid]

            for i in tip_id[str(cut_length)]:
                v = g_cut['tot'].parent(i)
                g_cut[str(cut_length)].property('k')[v] = g_cut[str(cut_length)].property('K')[v]

            g_cut[str(cut_length)] = flux.flux(g_cut[str(cut_length)], Jv = _Jv[count], psi_e = psi_base + DP_cnf[count], psi_base = psi_base,
                           invert_model = True, cut_and_flow = True)
            Jv = g_cut[str(cut_length)].property('J_out')[1]
            F += (Jv - _Jv[count]) ** 2.0

            count += 1


        return F

    def fun2(x):
        """
        Simulation of the flux at the different cut lengths according to the new parameter value

        Implementation 2: only axial_data is changed

        :param x: the array of adjusted parameters
        :return: F the sum((Jv - Jv_exp) ** 2.0)
        """
        # k0 = parameter.hydro['k0']

        axial_data = copy.deepcopy(parameter.hydro['axial_conductance_data'])
        axial_data[1] = list(x)

        g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], k_radial = k0 ,axial_data = axial_data, psi_base = psi_base,
                                                  psi_e = psi_base + DP_cnf[0])
        F = (Jv - _Jv[0])**2.0

        count = 1
        for cut_length in cut_n_flow_length:
            # _g = g_cut[str(cut_length)].copy() # not necessary and time-consuming

            for vid in g_cut[str(cut_length)].vertices_iter(g_cut['tot'].max_scale()):
                g_cut[str(cut_length)].property('K')[vid] = g_cut['tot'].property('K')[vid]
                g_cut[str(cut_length)].property('k')[vid] = g_cut['tot'].property('k')[vid]

            for i in tip_id[str(cut_length)]:
                v = g_cut['tot'].parent(i)
                g_cut[str(cut_length)].property('k')[v] = g_cut[str(cut_length)].property('K')[v]

            g_cut[str(cut_length)] = flux.flux(g_cut[str(cut_length)], Jv = _Jv[count], psi_e = psi_base + DP_cnf[count], psi_base = psi_base,
                           invert_model = True, cut_and_flow = True)
            Jv = g_cut[str(cut_length)].property('J_out')[1]
            F += (Jv - _Jv[count])**2.0

            count += 1

        return F

    def fun3(x):
        """
        Simulation of the flux at the different cut lengths according to the new parameter value

        Implementation 3: only k is changed

        :param x: the array of adjusted parameters
        :return: F the sum((Jv - Jv_exp) ** 2.0)
        """

        g_cut['tot']  = conductance.compute_k(g_cut['tot'] , k0 = x[0])

        g_cut['tot'] = flux.flux(g_cut['tot'], Jv =  _Jv[0], psi_e = psi_base + DP_cnf[0], psi_base = psi_base,
                           invert_model = True)
        Jv = g_cut['tot'].property('J_out')[1]
        F = (Jv - _Jv[0]) ** 2.0

        count = 1
        for cut_length in cut_n_flow_length:
            # _g = g_cut[str(cut_length)].copy() # not necessary and time-consuming

            for vid in g_cut[str(cut_length)].vertices_iter(g_cut['tot'].max_scale()):
                g_cut[str(cut_length)].property('K')[vid] = g_cut['tot'].property('K')[vid]
                g_cut[str(cut_length)].property('k')[vid] = g_cut['tot'].property('k')[vid]

            for i in tip_id[str(cut_length)]:
                v = g_cut['tot'].parent(i)
                g_cut[str(cut_length)].property('k')[v] = g_cut[str(cut_length)].property('K')[v]

            g_cut[str(cut_length)] = flux.flux(g_cut[str(cut_length)], Jv = _Jv[count], psi_e = psi_base + DP_cnf[count], psi_base = psi_base,
                               invert_model = True, cut_and_flow = True)
            Jv = g_cut[str(cut_length)].property('J_out')[1]
            F += (Jv - _Jv[count]) ** 2.0

            count += 1

        return F


    # architecture with filename in aqua team format
    if df_archi is None:
        # archi_f = glob.glob(parameter.archi['input_dir'] + parameter.archi['input_file'][0])
        # archi_f = archi_f[0]
        fname = str(Path(parameter.archi['input_dir']) / parameter.archi['input_file'][0]) # deals with path in windows
        archi_f = glob.glob(fname) # deals with wildcards ? deprecated ???
        archi_f = archi_f[0]

        df_archi = read_archi_data(archi_f) if parameter.archi['read_architecture'] else None
    index = archi_f.replace(glob.glob(parameter.archi['input_dir'])[0],"")

    # length law data: override if necessary
    if df_law is not None:
        parameter.archi['length_data'] = df_law

    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']

    columns = ['plant', 'cut length (m)', 'max_length', 'k (10-9 m/s/MPa)', 'length (m)', 'surface (m2)',
               'Jv (uL/s)', 'Jexp (uL/s)']

    results = {}
    for key in columns:
        results[key] = []


    # read the data measurements from data base
    if df_exp is not None:
        for key in df_exp['arch']:
            if str(key).lower() in index.lower():
                _list = df_exp[df_exp.arch == key].filter(regex = '^J').dropna(axis = 1).values.tolist()
                # parameter.exp['Jv'] = _list[0][0] # basal output flux full root (uncut)
                # _Jv = _list[0][1:]                # basal output flux cut root
                _Jv = _list[0]
                _list = df_exp[df_exp.arch == key].filter(regex = '^lcut').dropna(axis = 1).values.tolist()
                cut_n_flow_length = _list[0]      # cut lengthes
                _list = df_exp[df_exp.arch == key].filter(regex = '^dP').dropna(axis = 1).values.tolist()
                # the pressure difference is usually constant but sometimes, due to flow meter saturation, it may change
                # in that case a list of values is given
                if len(_list[0]) != 0:
                    DP_cnf = _list[0]
                else:
                    DP_cnf = []
                    DP_cnf.append(psi_e - psi_base) # for compatibility reason with first analysis on arabidopsis

                if len(DP_cnf) < len(cut_n_flow_length)+1: # if constant we create the list with the constant value
                    for i in range(1, len(cut_n_flow_length) + 1): DP_cnf.append(DP_cnf[0])
    else:
        _Jv = [parameter.exp['Jv']]
        cut_n_flow_length = []
        DP_cnf = [psi_e - psi_base]

    axfold = parameter.output['axfold'][0]
    radfold = parameter.output['radfold'][0]

    # g_cut['tot'], primary_length, _length, surface, seed = root_builder(df = df_archi, segment_length = parameter.archi['segment_length'],
    #     order_decrease_factor = parameter.archi['order_decrease_factor'], ref_radius = parameter.archi['ref_radius'])

    g_cut['tot'], primary_length, _length, surface, seed = \
        root_builder(df = df_archi,
                     primary_length = parameter.archi['primary_length'][0],
                    seed = parameter.archi['seed'][0],
                    delta = parameter.archi['branching_delay'][0],
                    nude_length = parameter.archi['nude_length'][0],
                    segment_length = parameter.archi['segment_length'],
                    length_data = parameter.archi['length_data'],
                    branching_variability = parameter.archi['branching_variability'],
                    order_max = parameter.archi['order_max'],
                    order_decrease_factor = parameter.archi['order_decrease_factor'],
                    ref_radius = parameter.archi['ref_radius'],
                    Flag_radius = Flag_radius)

    g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], psi_base = psi_base, psi_e = psi_base + DP_cnf[0])

    ###############################################################
    #### WARNING : the mtg property 'position' must stay unchanged
    ####           because the axial conductivity is placed according to it
    ###############################################################

    for cut_length in cut_n_flow_length:
        tip_id[str(cut_length)] = \
            flux.segments_at_length(g_cut['tot'], cut_length, dl = parameter.archi['segment_length'])
        g_cut[str(cut_length)] = \
            flux.cut_and_set_conductance(g_cut['tot'], cut_length, parameter.archi['segment_length'])
        # g_cut[str(cut_length)], surface = radius.compute_surface(g_cut[str(cut_length)])

    axial_data = list(conductance.axial(parameter.hydro['axial_conductance_data'], axfold))


    ###############################################################################################
    ## First adjustment: axfold, arfold that are coefficient factor of the radial conductivity k and
    ## and axial conductance K
    ###############################################################################################

    if Flag_Optim_K:
        # pre-optimization needed when there are several parameters so when K is optimized
        if Flag_verbose: print("*********** pre_optimization ************************")
        x = []
        bnds = []
        if Flag_Optim_K:
            x.append(axfold)
            bnds.append((1.0e-20, np.inf))
        if Flag_Optim_k:
            x.append(radfold)
            bnds.append((1.0e-20, np.inf))

        res = optimize.minimize(fun1, x, bounds = bnds, options = {'ftol': _tol})
        if Flag_Optim_k:
            radfold = res.x[-1] # always the last one even if the only one
            if Flag_verbose: print('pre-optimization ar: {:0.2e}'.format(res.x[-1]))
        if Flag_Optim_K:
            axfold = res.x[0]
            if Flag_verbose: print('pre-optimization ax: {:0.2e}'.format(res.x[0]))


            if Flag_verbose: print("****************************************************************")

    ## update the conductivities according to the first adjustment
    axial_data = list(conductance.axial(parameter.hydro['axial_conductance_data'], axfold))
    k0 = parameter.hydro['k0'] *radfold

    ###############################################################################################
    ## 2d adjustment:
    ##      -1 axial data adjusted
    ##      -2 radial conductivit adjusted
    ##      - 1 and 2 repeated until the k0 variation is below 0.1
    ###############################################################################################

    x = []
    x = axial_data[1]

    bnds = []
    n = len(x)
    for i, val in enumerate(x):
        bnds.append((1.0e-20, 1.0))
    # linear constraints lb <= A.dot(x) <= ub
    A = np.zeros((n, n))
    lb = np.full(n, -np.inf)
    ub = np.full(n, np.inf)
    l = parameter.hydro['axial_conductance_data'][0]
    a = dK_constraint # constraint on the 1st derivative
    ni = n - 1
    for i in range(ni): # downward derivative
        A[i, i] = -1.
        A[i, i + 1] = 1.
        lb[i] = a * (l[i+1]-l[i])
    i = ni
    A[i, i-1] = -1.
    A[i, i] = 1.
    lb[i] = a * (l[i]-l[i-1])


    k0_old = k0
    F_old = (Jv - _Jv[0])**2.0
    eps = 1e-9
    F = 1.
    if not (Flag_Optim_K or Flag_Optim_k):
        k0_old2 = k0
    else:
        k0_old2 = k0 + 10
    count2 = 0
    while abs(k0-k0_old2) > dk_max:
        k0_old2 = k0
        # parameter.hydro['k0'] = k0

# #ajout Mistral
#         _J_Sim=[]
#         lcut = []
#         lcut.append(0)
#
#         g_cut['tot']  = conductance.compute_k(g_cut['tot'] , k0 = k0)
#
#         g_cut['tot'] = flux.flux(g_cut['tot'], Jv =  _Jv[0], psi_e = psi_base + DP_cnf[0], psi_base = psi_base,
#                            invert_model = True)
#         _J_Sim.append(g_cut['tot'].property('J_out')[1])
#
#         count = 1
#         for cut_length in cut_n_flow_length:
#             # _g = g_cut[str(cut_length)].copy() # not necessary and time-consuming
#             lcut.append(g_cut['tot'].property('position')[1] - cut_length)
#             for vid in g_cut[str(cut_length)].vertices_iter(g_cut['tot'].max_scale()):
#                 g_cut[str(cut_length)].property('K')[vid] = g_cut['tot'].property('K')[vid]
#                 g_cut[str(cut_length)].property('k')[vid] = g_cut['tot'].property('k')[vid]
#
#             for i in tip_id[str(cut_length)]:
#                 v = g_cut['tot'].parent(i)
#                 g_cut[str(cut_length)].property('k')[v] = g_cut[str(cut_length)].property('K')[v]
#
#             g_cut[str(cut_length)] = flux.flux(g_cut[str(cut_length)], Jv = _Jv[count], psi_e = psi_base + DP_cnf[count], psi_base = psi_base,
#                                invert_model = True, cut_and_flow = True)
#             _J_Sim.append(g_cut[str(cut_length)].property('J_out')[1])
#
#
#             count += 1
#         plt.figure()
#         plt.scatter(lcut,_Jv, c = 'black')
#         plt.plot(lcut,_J_Sim, c = 'purple')
#         plt.savefig(str(count2) + '.png')
#         plt.close()
#         count2 += 1
# # end ajout Mistral
        ## -1 axial data adjusted
        #########################
        constraints = optimize.LinearConstraint(A, lb, ub)
        if Flag_Optim_K:
            res = optimize.minimize(fun2, x, bounds = bnds, constraints = constraints, options={'ftol': _tol})

            dKx = sum((x-res.x)**2.0)
            axial_data[1] = list(res.x)
            x = copy.deepcopy(res.x)

            if Flag_verbose: print('finished minimize K: [', ', '.join('{:0.2e}'.format(i) for i in res.x), ']')

        if Flag_Optim_k:
            ## -1 radial k adjusted
            #######################
            resk0 = optimize.minimize(fun3, k0, method = 'Nelder-Mead')


            if Flag_verbose: print('finished minimize k0: , {:0.2e}'.format(resk0.x[0])) #,
                                   # 'dk0**2.0 = {:0.2e}'.format((k0-resk0.x[0])**2.), 'dKx**2.0 = {:0.2e}'.format(dKx))

            k0 = resk0.x[0]
        else:
            k0_old2 = k0
    # print(resk0)
    ######################################
    ## Simulations with Kx and k adjusted
    ######################################

    primary_length = g_cut['tot'].property('position')[1]

    g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], k_radial = k0 ,axial_data = axial_data, psi_base = psi_base,
                                              psi_e = psi_base + DP_cnf[0])

    # add some properties for display
    g_cut['tot'].add_property('jsurf')
    length = g_cut['tot'].property('length')
    _radius = g_cut['tot'].property('radius')
    j = g_cut['tot'].property('j')
    jsurf = g_cut['tot'].property('jsurf')
    for v in g_cut['tot'].vertices_iter(scale = 1):
        jsurf[v] = j[v] / (_radius[v] * 2 * np.pi * length[v])
        
    results['plant'].append(index)
    results['max_length'].append(primary_length)
    results['cut length (m)'].append(0.0)
    results['k (10-9 m/s/MPa)'].append(k0)
    results['length (m)'].append(_length)
    results['surface (m2)'].append(surface)
    results['Jv (uL/s)'].append(Jv)
    results['Jexp (uL/s)'].append(_Jv[0])

    count = 1
    for cut_length in cut_n_flow_length:
        g_cut[str(cut_length)] = g_cut[str(cut_length)].copy()

        for vid in g_cut[str(cut_length)].vertices_iter(g_cut['tot'].max_scale()):
            g_cut[str(cut_length)].property('K')[vid] = g_cut['tot'].property('K')[vid]
            g_cut[str(cut_length)].property('k')[vid] = g_cut['tot'].property('k')[vid]

        for i in tip_id[str(cut_length)]:
            v = g_cut['tot'].parent(i)
            g_cut[str(cut_length)].property('k')[v] = g_cut[str(cut_length)].property('K')[v]

        g_cut[str(cut_length)] = flux.flux(g_cut[str(cut_length)], psi_e = psi_base + DP_cnf[count], psi_base = psi_base, invert_model = True)

        Jv = g_cut[str(cut_length)].property('J_out')[1]
        g_cut[str(cut_length)], surface = radius.compute_surface(g_cut[str(cut_length)])
        _length = g_cut[str(cut_length)].nb_vertices(scale = 1) * parameter.archi['segment_length']

        # add some properties for display
        g_cut[str(cut_length)].add_property('jsurf')
        length = g_cut[str(cut_length)].property('length')
        _radius = g_cut[str(cut_length)].property('radius')
        j = g_cut[str(cut_length)].property('j')
        jsurf = g_cut[str(cut_length)].property('jsurf')
        for v in g_cut[str(cut_length)].vertices_iter(scale = 1):
            jsurf[v] = j[v] / (_radius[v] * 2 * np.pi * length[v])

        primary_length = cut_length
        results['plant'].append(index)
        results['max_length'].append(primary_length)
        results['cut length (m)'].append(g_cut['tot'].property('position')[1] - primary_length)
        results['k (10-9 m/s/MPa)'].append(k0)
        results['length (m)'].append(_length)
        results['surface (m2)'].append(surface)
        results['Jv (uL/s)'].append(Jv)
        results['Jexp (uL/s)'].append(_Jv[count])
        count += 1


    dresults = pd.DataFrame(results, columns = columns)

    optim_results  = {}
    optim_results['x'] = copy.deepcopy(parameter.hydro['axial_conductance_data'][0])
    optim_results['K 1st'] = copy.deepcopy(parameter.hydro['axial_conductance_data'][1])

    if Flag_Optim_K:
        _x = list(res.x)
        optim_results['K optimized'] = copy.deepcopy(_x)
    else:
        optim_results['K optimized'] = optim_results['K 1st']

    doptim = pd.DataFrame(optim_results, columns = ['x', 'K 1st', 'K optimized'])

    df = pd.concat([dresults, doptim], axis = 1)

    if Flag_verbose:
        pd.set_option('display.max_columns', None)
        pd.set_option('display.expand_frame_repr', False)
        print(dresults)
        if Flag_Optim_K: print(doptim)

    if output is not None: df.to_csv(output, index = False)

    g_cut[0] = g_cut.pop('tot')
    icut = 1
    for cut_length in cut_n_flow_length:
        g_cut[icut] = g_cut.pop(str(cut_length))
        icut += 1

    return df, g_cut
