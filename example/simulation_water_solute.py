"""
Perform direct simulations using the water and solute transport module of Hydroroot.
May adjust k to fit the flux value given in the yaml file
Usage:
    %run adjustment_K_k_Js_Ps.py [-h] [-o OUTPUTFILE] [-op] [-v] [-d DATA] inputfile

        positional arguments:
            inputfile             yaml input file
        optional arguments:
            -h, --help            show this help message and exit
            -o , --outputfile     output csv file (default: out.csv)
            -op, --optimize       parameters adjustment (default: False)
            -v, --verbose         display parameter values during adjustment (default: False)

Inputs:
    - yaml file given in command line argument

Outputs:
    - console:
        - 'plant or seed', 'Jv (uL/s)', 'total length (m)', 'surface (m2)'
    - outputfile (csv):
        - columns = ['plant or seed', 'primary_length (m)', 'k (10-9 m/s/MPa)', 'total length (m)', 'surface (m2)', 'Jv (uL/s)']

"""
import glob
import argparse
import math
import time
import pandas as pd

from hydroroot.read_file import read_archi_data
from hydroroot.main import root_builder
from hydroroot import flux
from hydroroot.init_parameter import Parameters
from hydroroot.conductance import set_conductances, axial
from hydroroot.water_solute_transport import pressure_calculation_no_non_permeating_solutes, \
    init_some_MTG_properties, osmotic_p_peg

start_time = time.time()

###############################################################################
# Command line arguments
###############################################################################

parameter = Parameters()

parser = argparse.ArgumentParser(description='run direct simulation of water-solute HydroRoot, or adjust parameters on Jv(P) or Cut and flow or both data.')
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("-o", "--outputfile", default = 'out.csv', help="output csv file (default: out.csv)")
parser.add_argument("-op", "--optimize", default = False, help="parameters adjustment (default: False)", action="store_true")
parser.add_argument("-v", "--verbose", default = False, help="display parameter values during adjustment (default: False)", action="store_true")
args = parser.parse_args()
filename = args.inputfile
output = args.outputfile
Flag_Optim = args.optimize
Flag_verbose = args.verbose
parameter.read_file(filename)

###############################################################################
# Some Global boolean Flags: Flags I don't change often could be passed by command line
###############################################################################
Flag_radius = True  # True if radii furnished in architecture file used them, otherwise use ref_radius and so on from yaml file

###############################################################################
# Functions
###############################################################################

def Jv_calculation_non_permeating(g, axial_pr, k0_pr, Pe, Pbase, J_s = 0.0, sigma = 1.0, Ce = 0.0, P_s = 0.0, Cse = 0.0):
    """
    :Parameters:
    	- g: (MTG)
    	- axial_pr: (float list) axial conductance, K(x) data, 2 lists of float (microL.m/(s.Mpa))
    	- k: (float) radial conductivity (microL/(s.MPa.m**2))
    	- Pe: (float) external pressure (Mpa)
    	- Pbase: (float) base pressure (Mpa)
    	- J_s: (float) pumping rate (mol/(m2.s))
    	- P_s: (float) permeability coefficient  (m/s)
    	- sigma: (float) reflection coefficient
    	- Ce: (float) concentration of non-permeating solutes (mol/microL)
    	- Cse: (float) concentration of permeating solutes (mol/microL)
    :Returns:
    """
    g = set_conductances(g, axial_pr = axial_data, k0_pr = k0_pr)  # initialization
    g = flux.flux(g, psi_e = Pe, psi_base = Pbase)  # initialization
    g = init_some_MTG_properties(g, tau = J_s, Cini = Cse, t = 1)  # initialization
    nb_v = g.nb_vertices()
    Fdx = 1.0
    Fdx_old = 1.
    while Fdx > eps:
        g, dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g, sigma = Sigma, tau = J_s, Ce = Ce,
                                                                               Ps = P_s, Cse = Cse, Pe = Pe, Pbase = Pbase)
        Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
        if abs(Fdx - Fdx_old) < eps: break
        Fdx_old = Fdx
    Jv = g.property('J_out')[1]
    return Jv

###############################################################################
# Main
###############################################################################

if __name__ == '__main__':
    # files names of the architecture if reconstructed from a file
    # if not we just give a dummy name for the loop used to launch run
    filename = []
    parameter.archi['seed'] = [1]
    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))

    ############################
    # get value from yaml file
    ############################
    primary_length = parameter.archi['primary_length'][0]
    seed = parameter.archi['seed'][0]
    axfold = parameter.output['axfold'][0]
    radfold = parameter.output['radfold'][0]
    psi_base = parameter.exp['psi_base']
    psi_e = parameter.exp['psi_e']
    J_s = parameter.solute['J_s']
    P_s = parameter.solute['P_s']
    Cse = parameter.solute['Cse'] * 1e-9 # mol/m3 -> mol/microL, external permeating solute concentration
    Ce = parameter.solute['Ce'] * 1e-9 # mol/m3 -> mol/microL, external non-permeating solute concentration
    Cini = Cse # initialization solute concentration into the xylem vessels
    Cpeg_ini = Ce # initialization non-permeating solute concentration into the xylem vessels: not 0.0 because more num instability
    Sigma = parameter.solute['Sigma'] # reflection coefficient, fixed in this script
    Pi_e_peg = osmotic_p_peg(Ce, unit_factor = 8.0e6)  # from Ce mol/microL to g/g, external osmotic pressure of non-permeating in MPa

    # Conductancies: mananging the fact there are or not different values between the primary and laterals
    # and the fact there are multiply by axfold and radfold
    k = parameter.hydro['k0'] * radfold
    axial_data = parameter.hydro['axial_conductance_data']
    axial_data = list(axial(axial_data, axfold))

    data = None
    row = None
    col = None

    eps = 1.0e-9 # global: stop criterion for the Newton-Raphson loop in Jv_P_calculation and Jv_cnf_calculation


    columns = ['plant or seed', 'primary_length (m)', 'k (10-9 m/s/MPa)', 'total length (m)', 'surface (m2)', 'Jv (uL/s)']
    results = {}
    for key in columns:
        results[key] = []

    # we can give a list of architecture file names
    for f in filename:
        # architecture file to dataframe
        df = read_archi_data(f) if parameter.archi['read_architecture'] else None

        # build the MTG
        ################
        g, primary_length, total_length, surface, _seed = root_builder(primary_length, df = df, seed = seed,
                                           segment_length = parameter.archi['segment_length'],
                                           length_data = parameter.archi['length_data'],
                                           branching_variability = parameter.archi['branching_variability'],
                                           order_max = parameter.archi['order_max'],
                                           order_decrease_factor = parameter.archi['order_decrease_factor'],
                                           ref_radius = parameter.archi['ref_radius'],
                                           Flag_radius = Flag_radius)

        # run simulation or optimization
        #################################
        Jv = Jv_calculation_non_permeating(g, axial_pr = axial_data, k0_pr = k, Pe = psi_e, Pbase = psi_base,
                                                J_s = J_s, sigma = Sigma, Ce = Ce, P_s = P_s, Cse = Cse)

        # simple Newton-Raphson loop to adjust k
        if Flag_Optim:
            k_old = k
            F_old = (Jv - parameter.exp['Jv'])**2.0
            k *= 0.9
            eps = 1e-15
            F = 1.
            # Newton-Raphson loop to get k
            while (F > eps):
                Jv = Jv_calculation_non_permeating(g, axial_pr = axial_data, k0_pr = k, Pe = psi_e,
                                                        Pbase = psi_base,
                                                        J_s = J_s, sigma = Sigma, Ce = Ce, P_s = P_s, Cse = Cse)

                F = (Jv - parameter.exp['Jv']) ** 2.0

                if Flag_verbose: print(k)

                if abs(F) > eps:
                    dfdk = (F - F_old) / (k - k_old)
                    k_old = k
                    k = k_old - F / dfdk
                    while k < 1.0e-3:
                        k = 0.5 * k_old
                    F_old = F

        if df is None:
            index = _seed
        else:
            index = f.replace(glob.glob(parameter.archi['input_dir'])[0], "")
        results['plant or seed'].append(index)
        results['primary_length (m)'].append(primary_length)
        results['k (10-9 m/s/MPa)'].append(k)
        results['total length (m)'].append(total_length)
        results['surface (m2)'].append(surface)
        results['Jv (uL/s)'].append(Jv)

    dr = pd.DataFrame(results, columns = columns)
    print(dr.loc[:, ['plant or seed', 'Jv (uL/s)', 'total length (m)', 'surface (m2)']])
    if output is not None:
        dr.to_csv(output,  index = False)
    print('running time is ', time.time() - start_time)
