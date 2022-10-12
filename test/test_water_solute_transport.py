import math

import sys
sys.path.insert(0, '../src')

from hydroroot.conductance import set_conductances
from hydroroot.main import root_builder, flux
from hydroroot.init_parameter import Parameters
from hydroroot.water_solute_transport import viscosity_peg, osmotic_p_peg, derivative_viscosity_peg, \
    derivative_osmotic_p_peg, pressure_calculation, init_some_MTG_properties, pressure_calculation_no_non_permeating_solutes

def test_viscosity_etc():
    # test function that calculate viscosity, osmotic pressure and their derivatives according to non-permeating concentrations
    eps = 1.0e-15
    C = 1.0
    ufactor = 1.0e-2
    w = C * ufactor
    T = 25.0
    mu = viscosity_peg(C, unit_factor = ufactor)
    mu2 = -17.4 + 18.4 * math.exp(w/0.279)

    assert abs(mu - mu2) < eps

    p = osmotic_p_peg(C, unit_factor = ufactor)
    p2 = (1.29 * w ** 2 * T - 140 * w ** 2 - 4.0 * w) * 0.1

    assert abs(p - p2) < eps

    dmu = derivative_viscosity_peg(C, unit_factor = ufactor)
    dmu2 = 18.4 * math.exp(w/0.279) / 0.279

    assert abs(dmu - dmu2) < eps

    dp = derivative_osmotic_p_peg(C, unit_factor = ufactor)
    dp2 = (2.58 * w  * T - 280 * w  - 4.0 ) * 0.1 * ufactor

    assert abs(dp - dp2) < eps

def test_pressure_calculation():
    # test pressure_calculation and pressure_calculation_no_non_permeating_solutes: compare resulting flux when Ce=0
    #   Ce: non-permeating solute concentration
    # test also root_builder and init_some_MTG_properties
    parameter = Parameters()
    parameter.read_file('data/parameters_test_yaml.yml')

    g, _p, _l, _s, _seed = root_builder(0.1, segment_length = 1.0e-3,
                                        length_data = parameter.archi['length_data'],
                                        order_max = 0,
                                        order_decrease_factor = parameter.archi['order_decrease_factor'],
                                        ref_radius = parameter.archi['ref_radius'])

    # needed pressure_calculation
    g.add_property('C')
    g.add_property('Cpeg')
    g.add_property('JC')
    g.add_property('pi')
    g.add_property('theta')
    g.add_property('DPsi')
    g.add_property('J_s')
    g.add_property('js')

    J_s = parameter.solute['J_s']
    P_s = parameter.solute['P_s']
    Cini = Cse = parameter.solute['Cse']
    Ce = parameter.solute['Ce']
    Sigma = parameter.solute['Sigma']
    psi_base = parameter.exp['psi_base']
    psi_e = parameter.exp['psi_e']
    eps = 1.0e-9

    g = set_conductances(g, axial_pr = parameter.hydro['axial_conductance_data'], k0_pr = parameter.hydro['k0'])
    g = flux.flux(g, psi_e = psi_e, psi_base = psi_base)
    g = init_some_MTG_properties(g, tau = J_s, Cini = Cini, t = 1, Ps = P_s)
    nb_v = g.nb_vertices()
    Fdx = 1.0
    Fdx_old = 1.
    dP = psi_e - psi_base
    while Fdx > eps:
        g, dx, data, row, col = pressure_calculation(g, sigma = Sigma, Ce = Ce, Cse = Cse, dP = dP)
        Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
        local_j = g.property('J_out')[1]
        if abs(Fdx - Fdx_old) < eps: break
        Fdx_old = Fdx

    g = set_conductances(g, axial_pr = parameter.hydro['axial_conductance_data'], k0_pr = parameter.hydro['k0'])
    g = flux.flux(g, psi_e = psi_e, psi_base = psi_base)
    g = init_some_MTG_properties(g, tau = J_s, Cini = Cini, t = 1)
    nb_v = g.nb_vertices()
    Fdx = 1.0
    Fdx_old = 1.
    dP = psi_e - psi_base
    while Fdx > eps:
        g, dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g, sigma = Sigma, Ce = Ce, Cse = Cse, dP = dP)
        Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
        local_j2 = g.property('J_out')[1]
        if abs(Fdx - Fdx_old) < eps: break
        Fdx_old = Fdx

    assert abs(local_j-local_j2) < eps