###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to adjust axial conductance data and radial conductivity on
#   cut and flow experiments
###############################################################################

######
# Imports

# VERSION = 2

import numpy as np
import pandas as pd
import glob
import copy
import argparse

from openalea.mtg import traversal
from scipy import optimize

from hydroroot import radius, flux, conductance
from hydroroot.generator.measured_root import mtg_from_aqua_data
from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters

import datetime, time

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

# read architecture file
def read_archi_data(fn):
    """
    Read a csv (tab separated) file with the architecture in the following format
        |'distance_from_base_(mm)' | 'lateral_root_length_(mm)' | order |
        |float | float | string|
        order = 1 for laterals of 1st order ob the primary
        order = n-m for the lateral number m on the lateral number n of 1st order
        order = n-m-o for the lateral number o of the previous one
        etc.
        Each branch finish with a nude part, i.e. a distance from base (the tip value) and a zero length

    :param fn: string - the architecture filename in csv format

    :return: DataFrame
    """
    df = pd.read_csv(fn, sep = '\t', dtype = {'order': str})
    df['db'] = df['distance_from_base_(mm)'] * 1.e-3
    df['lr'] = df['lateral_root_length_(mm)'] * 1.e-3
    return df

def radial(v = 92, acol = [], scale = 1):
    xr = acol[0]  # at this stage kr constant so the same x than Ka
    yr = [v * scale] * len(xr)
    return xr, yr

def axial(acol = [], scale = 1):
    x, y = acol
    y = [a * scale for a in y]
    return x, y

def root_creation(df = None):
    """
    creation of an mtg with properties like radius and vertex length set.

    The MTG is either generated or created from a data.
    The radius and vertex length properties are set.
    The following properties are computed: length, position, mylength, surface, volume, total length,
        primary root length, nb of intercepts

    :param:
        primary_length: primary root length for generated mtg
        seed:  seed for generated mtg, if None randomly generated
        delta: branching delay  for generated mtg
        nude_length: length from tip without lateral for generated mtg
        df: pandas DataFrame with the architecture data to be reconstructed
    :return:
        g: MTG with the different properties set or computed (see comments above),
        primary_length: primary root length (output for generated mtg)
        _length: total root length
        surface: total root surface
    """
    g = mtg_from_aqua_data(df, parameter.archi['segment_length'])


    # compute radius property on MTG
    g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])

    # compute length property and parametrisation
    g = radius.compute_length(g, parameter.archi['segment_length'])
    g = radius.compute_relative_position(g)

    # Calculation of the distance from base of each vertex, used for cut and flow
    # Remark: this calculation is done in flux.segments_at_length; analysis.nb_roots but there is a concern with the
    # parameter dl which should be equal to vertex length but which is not pass
    _mylength = {}
    for v in traversal.pre_order2(g, 1):
        pid = g.parent(v)
        _mylength[v] = _mylength[pid] + parameter.archi['segment_length'] if pid else parameter.archi['segment_length']
    g.properties()['mylength'] = _mylength

    # _length is the total length of the RSA (sum of the length of all the segments)
    _length = g.nb_vertices(scale = 1) * parameter.archi['segment_length']
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    v_base = g.component_roots_at_scale_iter(g.root, scale = g.max_scale()).next()
    primary_length = g.property('position')[v_base]


    return g, primary_length, _length, surface

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None):
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

def fun1(x):
    """
    Simulation of the flux at the different cut lengths according to the new parameter value

    Implementation 1: only axfold (Kx factor) and radfold (k radial factor) are changed

    :param x: the array of adjusted parameters
    :return: F the sum((Jv - Jv_exp) ** 2.0)
    """
    axfold = x[0]
    radfold = x[-1]

    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']

    g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], radfold = radfold, axfold = axfold)
    F = (Jv - parameter.exp['Jv']) ** 2.0
    count = 0
    for cut_length in cut_n_flow_length:
        _g = g_cut[str(cut_length)].copy()

        for vid in _g.vertices_iter(g_cut['tot'].max_scale()):
            _g.property('K')[vid] = g_cut['tot'].property('K')[vid]
            _g.property('k')[vid] = g_cut['tot'].property('k')[vid]

        for i in tip_id[str(cut_length)]:
            v = g_cut['tot'].parent(i)
            _g.property('k')[v] = _g.property('K')[v]

        _g = flux.flux(_g, Jv = _Jv[count], psi_e = psi_e, psi_base = psi_base,
                       invert_model = True, cut_and_flow = True)
        keq = _g.property('Keq')[1]
        Jv = keq * (psi_e - psi_base)
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
    k0 = parameter.hydro['k0']

    axial_data = copy.deepcopy(parameter.hydro['axial_conductance_data'])
    axial_data[1] = list(x)

    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']

    g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], k_radial = k0 ,axial_data = axial_data)
    F = (Jv - parameter.exp['Jv'])**2.0 
    
    count = 0
    for cut_length in cut_n_flow_length:
        _g = g_cut[str(cut_length)].copy()

        for vid in _g.vertices_iter(g_cut['tot'].max_scale()):
            _g.property('K')[vid] = g_cut['tot'].property('K')[vid]
            _g.property('k')[vid] = g_cut['tot'].property('k')[vid]

        for i in tip_id[str(cut_length)]:
            v = g_cut['tot'].parent(i)
            _g.property('k')[v] = _g.property('K')[v]

        _g = flux.flux(_g, Jv = _Jv[count], psi_e = psi_e, psi_base = psi_base,
                       invert_model = True, cut_and_flow = True)
        keq = _g.property('Keq')[1]
        Jv = keq * (psi_e - psi_base)
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
    count = 0
    g_cut['tot']  = conductance.compute_k(g_cut['tot'] , k0 = x[0])
    
    g_cut['tot'] = flux.flux(g_cut['tot'], Jv =  parameter.exp['Jv'], psi_e = psi_e, psi_base = psi_base,
                       invert_model = True)
    keq = g_cut['tot'].property('Keq')[1]
    Jv = keq * (psi_e - psi_base)
    F = (Jv - parameter.exp['Jv']) ** 2.0

    for cut_length in cut_n_flow_length:
        _g = g_cut[str(cut_length)].copy()

        for vid in _g.vertices_iter(g_cut['tot'].max_scale()):
            _g.property('K')[vid] = g_cut['tot'].property('K')[vid]
            _g.property('k')[vid] = g_cut['tot'].property('k')[vid]

        for i in tip_id[str(cut_length)]:
            v = g_cut['tot'].parent(i)
            _g.property('k')[v] = _g.property('K')[v]

        _g = flux.flux(_g, Jv = _Jv[count], psi_e = psi_e, psi_base = psi_base,
                       invert_model = True, cut_and_flow = True)
        keq = _g.property('Keq')[1]
        Jv = keq * (psi_e - psi_base)
        F += (Jv - _Jv[count]) ** 2.0

        count += 1

    return F

if __name__ == '__main__':


    dK_constraint = -3.e-2 # dK/dx >= dK_constraint
    _tol = 1.0e-9
    
    filename = []
    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))

    fn = 'data/arabido_data.csv'
    df_exp = pd.read_csv(fn, sep = ',', keep_default_na = True)

    # predict the number of simulation run
    nb_steps = len(filename)
    print 'Simulation runs: ', nb_steps
    print '#############################'


    psi_e = parameter.exp['psi_e']
    psi_base = parameter.exp['psi_base']

    columns = ['plant', 'cut length (m)', 'primary_length (m)', 'k (10-8 m/s/MPa)', '_length (m)', 'surface (m2)', 'Jv (uL/s)', 'Jexp (uL/s)']

    results = {}
    for key in columns:
        results[key] = []

    for f in filename:
        df = read_archi_data(f) if parameter.archi['read_architecture'] else None
        index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")

        # read the data measurements from data base
        for key in df_exp['arch']:
            if key in index.lower():
                _list = df_exp[df_exp.arch == key].filter(regex = '^J').dropna(axis = 1).values.tolist()
                parameter.exp['Jv'] = _list[0][0] # basal output flux full root (uncut)
                _Jv = _list[0][1:]                # basal output flux cut root
                _list = df_exp[df_exp.arch == key].filter(regex = '^lcut').dropna(axis = 1).values.tolist()
                cut_n_flow_length = _list[0]      # cut lengthes

        axfold = parameter.output['axfold'][0]
        radfold = parameter.output['radfold'][0]

        g_cut['tot'], primary_length, _length, surface = root_creation(df = df)

        g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'])

        ###############################################################
        #### WARNING : the mtg property 'position' must stay unchanged
        ####           because the axial conductivity is placed according to it
        ###############################################################

        for cut_length in cut_n_flow_length:
            tip_id[str(cut_length)] = \
                flux.segments_at_length(g_cut['tot'], cut_length, dl = parameter.archi['segment_length'])
            g_cut[str(cut_length)] = \
                flux.cut_and_set_conductance(g_cut['tot'], cut_length, parameter.archi['segment_length'])
            g_cut[str(cut_length)], surface = radius.compute_surface(g_cut[str(cut_length)])

        axial_data = list(axial(parameter.hydro['axial_conductance_data'], axfold))


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

            print "finished minimize ax, ar", res
            print "*******************************************************************************"
        
        ## update the conductivities according to the first adjustment
        axial_data = list(axial(parameter.hydro['axial_conductance_data'], axfold))
        parameter.hydro['k0'] = parameter.hydro['k0'] *radfold
        
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

        k0 = parameter.hydro['k0']
        k0_old = k0                        # count = 0

        F_old = (Jv - parameter.exp['Jv'])**2.0

        eps = 1e-9
        F = 1.
        if not Flag_Optim:
            k0_old2 = k0
        else:
            k0_old2 = k0 + 10
        while abs(k0-k0_old2) > 1.0e-1:
            k0_old2 = k0
            parameter.hydro['k0'] = k0

            ## -1 axial data adjusted
            #########################
            constraints = optimize.LinearConstraint(A, lb, ub)
            res = optimize.minimize(fun2, x, bounds = bnds, constraints = constraints, options={'ftol': _tol})

            dKx = sum((x-res.x)**2.0)
            axial_data[1] = list(res.x)
            x = copy.deepcopy(res.x)

            print "finished minimize Kx", res

            ## -1 radial k adjusted
            #######################
            resk0 = optimize.minimize(fun3, k0, method = 'Nelder-Mead')

            print 'Simu, ', k0, resk0.fun, resk0.x[0], 'dk0 = ', (k0-resk0.x[0])**2., 'dKx = ', dKx

            k0 = resk0.x[0]
        

        parameter.hydro['k0'] = k0

    primary_length = g_cut['tot'].property('position')[1]
    g_cut['tot'], Keq, Jv = hydro_calculation(g_cut['tot'], k_radial = k0 ,axial_data = axial_data)

    results['plant'].append(index)
    results['primary_length (m)'].append(primary_length)
    results['cut length (m)'].append(0.0)
    results['k (10-8 m/s/MPa)'].append(k0*0.1) #uL/s/MPa/m2 -> 10-8 m/s/MPa
    results['_length (m)'].append(_length)
    results['surface (m2)'].append(surface)
    results['Jv (uL/s)'].append(Jv)
    results['Jexp (uL/s)'].append(parameter.exp['Jv'])
    
    print primary_length, Jv

    ######################################
    ## Simulations with Kx and k adjusted
    ######################################
    count = 0
    for cut_length in cut_n_flow_length:
        _g = g_cut[str(cut_length)].copy()

        for vid in _g.vertices_iter(g_cut['tot'].max_scale()):
            _g.property('K')[vid] = g_cut['tot'].property('K')[vid]
            _g.property('k')[vid] = g_cut['tot'].property('k')[vid]

        for i in tip_id[str(cut_length)]:
            v = g_cut['tot'].parent(i)
            _g.property('k')[v] = _g.property('K')[v]

        _g = flux.flux(_g, psi_e = psi_e, psi_base = psi_base, invert_model = True)

        Jv = _g.property('J_out')[1]
        _g, surface = radius.compute_surface(_g)
        _length = _g.nb_vertices(scale = 1) * parameter.archi['segment_length']

        primary_length = cut_length
        results['plant'].append(index)
        results['primary_length (m)'].append(primary_length)
        results['cut length (m)'].append(g_cut['tot'].property('position')[1] - primary_length)
        results['k (10-8 m/s/MPa)'].append(k0*0.1) #uL/s/MPa/m2 -> 10-8 m/s/MPa
        results['_length (m)'].append(_length)
        results['surface (m2)'].append(surface)
        results['Jv (uL/s)'].append(Jv)
        results['Jexp (uL/s)'].append(_Jv[count])
        count += 1
        print index, primary_length, k0*0.1, _length, surface, Jv
        

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
    else:
        df = dresults

    if output is not None: df.to_csv(output, index = False)
    print 'running time is ', time.time() - start_time