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
    - console:
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

import numpy as np
import glob
import argparse
import time
import copy
import pandas as pd
import matplotlib.pyplot as plt
# import k3d
from scipy import optimize

from openalea.plantgl.all import *
from openalea.mtg import *
from oawidgets.plantgl import PlantGL

from hydroroot import flux, conductance, radius
from hydroroot.main import hydroroot_flow, root_builder
from hydroroot.init_parameter import Parameters
from hydroroot.read_file import read_archi_data
from hydroroot.display import mtg_scene
from ipywidgets import interact, fixed

from openalea.plantgl.algo.view import view


results = {}
g = None
Jv_global = 1.0
g_cut = {}
tip_id = {}
cut_n_flow_length = []

start_time = time.time()

parameter = Parameters()

################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("-o", "--outputfile", help="output csv file")
parser.add_argument("-op", "--optimize", help="optimize K and k", action="store_true")
# parser.add_argument("-v", "--verbose", help="display convergence", action="store_true")
args = parser.parse_args()
filename = args.inputfile
output = args.outputfile
if output is None: output = "out.csv"
Flag_Optim = args.optimize
if Flag_Optim is None: Flag_Optim = False
# Flag_verbose = args.verbose
# if Flag_verbose is None: Flag_verbose = False
parameter.read_file(filename)


def plot_architecture():

    def my_view(cut: str = 'tot', prop: str = 'j', imgsize: tuple = (800, 800), perspective: bool = True,
                zoom: float = 1,
                azimuth: float = 0, elevation: float = 0, line_width = 1.0):
        list_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K', 'axial flux': 'J_out'}

        g = g_cut[cut].copy()

        keys = list(g.property('radius').keys())
        radius = np.array(list(g.property('radius').values()))
        new_radius = radius * line_width
        g.properties()['radius'] = dict(list(zip(keys, new_radius)))

        s = mtg_scene(g, prop_cmap = list_prop[prop], has_radius = True)
        return view(scene = s, imgsize = imgsize, perspective = perspective, zoom = zoom, azimuth = azimuth, elevation = elevation)

    _list = ['tot']
    for i in cut_n_flow_length:
        _list.append(str(i))

    interact(my_view, cut = _list, prop = ['radial flux', 'xylem pressure', 'axial K', 'axial flux'],
             imgsize = fixed((800, 800)),
             perspective = False, zoom = (0.01, 1), azimuth = (-180, 180), elevation = (-90, 90), line_width = (1, 5))


def plot_3D():
    def my_view(cut: str = 'tot', prop: str = 'j', line_width = 1.0):
        list_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K', 'axial flux': 'J_out'}

        g = g_cut[cut].copy()

        keys = list(g.property('radius').keys())
        radius = np.array(list(g.property('radius').values()))
        new_radius = radius * line_width
        g.properties()['radius'] = dict(list(zip(keys, new_radius)))

        s = mtg_scene(g, prop_cmap = list_prop[prop], has_radius = True)

        return PlantGL(s)

    _list = ['tot']
    for i in cut_n_flow_length:
        _list.append(str(i))

    interact(my_view, cut = _list, prop = ['radial flux', 'xylem pressure', 'axial K', 'axial flux'],
             line_width = (1, 5))

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None, psi_base = 0.1, psi_e = 0.1):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    Kexp_axial_data = conductance.axial(axial_data, axfold)
    k_radial_data = conductance.radial(k_radial, axial_data, radfold)
    # k_radial_data = conductance.radial_step(k_radial,3.0,x_step = 0.01, dx = parameter.archi['segment_length'], scale = radfold)

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
        # g_cut[str(cut_length)] = g_cut[str(cut_length)].copy()

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
        # g_cut[str(cut_length)] = g_cut[str(cut_length)].copy()

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
        # g_cut[str(cut_length)] = g_cut[str(cut_length)].copy()

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

if __name__ == '__main__':


    dK_constraint = -3.e-2 # dK/dx >= dK_constraint # arabidopsis paper Boursiac2022
    _tol = 1.0e-9
    
    filename = []
    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))

    fn = 'data/cnf_data.csv'
    df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)

    # predict the number of simulation run
    nb_steps = len(filename)
    print('Simulation runs: ', nb_steps)
    print('#############################')


    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']

    columns = ['plant', 'cut length (m)', 'primary_length (m)', 'k (10-9 m/s/MPa)', '_length (m)', 'surface (m2)', 'Jv (uL/s)', 'Jexp (uL/s)']

    results = {}
    for key in columns:
        results[key] = []

    for f in filename:
        df = read_archi_data(f) if parameter.archi['read_architecture'] else None
        index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")

        # read the data measurements from data base
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
                    if len(DP_cnf) < len(cut_n_flow_length)+1: # if constant we create the list with the constant value
                        for i in range(1, len(cut_n_flow_length) + 1): DP_cnf.append(_list[0][0])

        axfold = parameter.output['axfold'][0]
        radfold = parameter.output['radfold'][0]

        g_cut['tot'], primary_length, _length, surface, seed = root_builder(df = df, segment_length = parameter.archi['segment_length'],
            order_decrease_factor = parameter.archi['order_decrease_factor'], ref_radius = parameter.archi['ref_radius'])

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
        if Flag_Optim:
            optim_k0 = True
            optim_K = True
            x = []
            if optim_K: x.append(axfold)
            if optim_k0: x.append(radfold)

            bnds = []
            if optim_K: bnds.append((1.0e-20, np.inf))
            if optim_k0: bnds.append((1.0e-20, np.inf))
            res = optimize.minimize(fun1, x, bounds = bnds, options = {'ftol': _tol})
            radfold = res.x[-1] # always the last one even if the only one
            axfold = res.x[0]

            print("finished minimize ax, ar", res.x)
            print("**********************************")
        
        ## update the conductivities according to the first adjustment
        axial_data = list(conductance.axial(parameter.hydro['axial_conductance_data'], axfold))
        k0 = parameter.hydro['k0'] *radfold
        
        ###############################################################################################
        ## 2d adjustment: 
        ##      -1 axial data adjusted
        ##      -2 radial conductivit adjusted
        ##      - 1 and 2 repeated until the k0 variation is below 0.1
        ###############################################################################################
        
        optim_k0 = False
        optim_K = True
        x = []
        x = axial_data[1] #copy.deepcopy(parameter.hydro['axial_conductance_data'][1])

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
        if not optim_k0:
            ni = n - 1
        else:
            ni = n - 2
        for i in range(ni): # downward derivative
            A[i, i] = -1.
            A[i, i + 1] = 1.
            lb[i] = a * (l[i+1]-l[i])

        i = ni
        A[i, i-1] = -1.
        A[i, i] = 1.
        lb[i] = a * (l[i]-l[i-1])

        # k0 = parameter.hydro['k0']
        k0_old = k0                        # count = 0

        F_old = (Jv - _Jv[0])**2.0

        eps = 1e-9
        F = 1.
        if not Flag_Optim:
            k0_old2 = k0
        else:
            k0_old2 = k0 + 10
        while abs(k0-k0_old2) > 1.0e-1:
            k0_old2 = k0
            # parameter.hydro['k0'] = k0

            ## -1 axial data adjusted
            #########################
            constraints = optimize.LinearConstraint(A, lb, ub)
            res = optimize.minimize(fun2, x, bounds = bnds, constraints = constraints, options={'ftol': _tol})

            dKx = sum((x-res.x)**2.0)
            axial_data[1] = list(res.x)
            x = copy.deepcopy(res.x)

            # print("finished minimize Kx", res.x)

            ## -1 radial k adjusted
            #######################
            resk0 = optimize.minimize(fun3, k0, method = 'Nelder-Mead')

            print('Simu, ', k0, resk0.fun, resk0.x[0], 'dk0 = ', (k0-resk0.x[0])**2., 'dKx = ', dKx)

            k0 = resk0.x[0]
        
        print('**********************')
        print('optimization finished')
        print('**********************')

        # parameter.hydro['k0'] = k0

    primary_length = g_cut['tot'].property('position')[1]

    g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], k_radial = k0 ,axial_data = axial_data, psi_base = psi_base,
                                              psi_e = psi_base + DP_cnf[0])

    results['plant'].append(index)
    results['primary_length (m)'].append(primary_length)
    results['cut length (m)'].append(0.0)
    results['k (10-9 m/s/MPa)'].append(k0)
    results['_length (m)'].append(_length)
    results['surface (m2)'].append(surface)
    results['Jv (uL/s)'].append(Jv)
    results['Jexp (uL/s)'].append(_Jv[0])
    
    print(index, primary_length, k0, _length, surface, Jv)

    ######################################
    ## Simulations with Kx and k adjusted
    ######################################
    count = 1
    for cut_length in cut_n_flow_length:
        # g_cut[str(cut_length)] = g_cut[str(cut_length)].copy()

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

        primary_length = cut_length
        results['plant'].append(index)
        results['primary_length (m)'].append(primary_length)
        results['cut length (m)'].append(g_cut['tot'].property('position')[1] - primary_length)
        results['k (10-9 m/s/MPa)'].append(k0)
        results['_length (m)'].append(_length)
        results['surface (m2)'].append(surface)
        results['Jv (uL/s)'].append(Jv)
        results['Jexp (uL/s)'].append(_Jv[count])
        count += 1
        print(index, primary_length, k0, _length, surface, Jv)
        

    dresults = pd.DataFrame(results, columns = columns)

    ax = dresults.plot.scatter('cut length (m)', 'Jexp (uL/s)', c = 'black')
    dresults.plot.line('cut length (m)', 'Jv (uL/s)', c = 'purple', ax = ax)


    if Flag_Optim:
        optim_results  = {}
        optim_results['x'] = copy.deepcopy(parameter.hydro['axial_conductance_data'][0])
        optim_results['K 1st'] = copy.deepcopy(parameter.hydro['axial_conductance_data'][1])
        _x = list(res.x)
        optim_results['K optimized'] = copy.deepcopy(_x)

        doptim = pd.DataFrame(optim_results, columns = ['x', 'K 1st', 'K optimized'])
        df = pd.concat([dresults, doptim], axis = 1)

        ax_K = doptim.plot.line('x', 'K 1st', c = 'black')
        doptim.plot.line('x', 'K optimized', c = 'purple', ax = ax_K)
        
        print('************************************************************')
        print('radial k 1st: ', parameter.hydro['k0'], ', k optimized: ', k0)
        print('************************************************************')
        fig, ax2 = plt.subplots()
        f1=ax2.bar(x=['radial k 1st','k optimized'], height = [parameter.hydro['k0'],k0])
        ax2.set_ylabel('radial k (10-9 m/s/MPa)')
    else:
        df = dresults
    
    

    plt.show()

    # plot_architecture()
    
    if output is not None: df.to_csv(output, index = False)
    print('running time is ', time.time() - start_time)
