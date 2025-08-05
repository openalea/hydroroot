
from openalea.hydroroot import flux
from openalea.hydroroot.main import root_builder
from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.read_file import read_archi_data
from openalea.hydroroot.conductance import set_conductances
from openalea.mtg import traversal

def test_implementation():
    '''
    pure water solver i.e. electical network
    test if flux.flux gives the same results with has_soil=True when setting the same value to each psi_e[v] and
    has_soil=False (i.e. psi_e=None)
    test also if the equivalent calculation of K and Pext gives the same
    '''
    parameter = Parameters()
    parameter.read_file('data/parameters_test_yaml.yml')

    fname = parameter.archi['input_dir'] + parameter.archi['input_file'][0]
    df = read_archi_data(fname)
    g, primary_length, total_length, surface, seed = root_builder( primary_length = parameter.archi['primary_length'],
                                                                    delta = parameter.archi['branching_delay'],
                                                                    nude_length = parameter.archi['nude_length'],
                                                                    df = df,
                                                                    segment_length = parameter.archi['segment_length'],
                                                                    length_data = parameter.archi['length_data'],
                                                                    order_max = parameter.archi['order_max'],
                                                                    order_decrease_factor = parameter.archi['order_decrease_factor'],
                                                                    ref_radius = parameter.archi['ref_radius'])

    psi_e = g.property('psi_e')
    _h = g.property('base_length')
    hmax = max(_h.values())
    Pbase = parameter.exp['psi_base']
    Pe = parameter.exp['psi_e'] = 0.3
    # Select the base of the root
    v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))

    for v in traversal.post_order2(g, v_base):
        psi_e[v] = Pe #Pbase + (Pe - Pbase) / hmax * _h[v]

    g = set_conductances(g, axial_pr = parameter.hydro['axial_conductance_data'], k0_pr = parameter.hydro['k0'])
    g = flux.flux(g, psi_base = Pbase)  # flux.flux will use flux.has_soil = True

    psi_out = g.property('psi_out')
    psi_in = g.property('psi_in')
    Peq = g.property('Peq')
    Keq = g.property('Keq')
    Jv1 = - g.property('K')[v_base] * (psi_out[v_base] - psi_in[v_base])
    Jv2 = Keq[v_base] * (Peq[v_base] - Pbase)

    psi_e = None # flux.flux will use flux.has_soil = False
    g = flux.flux(g, psi_e = Pe, psi_base = Pbase)  # initialization
    Jv3 = - g.property('K')[v_base] * (g.property('psi_out')[v_base] - g.property('psi_in')[v_base])

    eps = 1.0e-10
    assert abs(Jv1-Jv2) <= eps
    assert abs(Jv1-Jv3) <= eps

    # now test if, with a gradient in the soil we have indeed a decreasing flux
    psi_e = g.property('psi_e')
    for v in traversal.post_order2(g, v_base):
        psi_e[v] = Pbase + (Pe - Pbase) / hmax * _h[v]
    g = flux.flux(g, psi_base = Pbase)  # flux.flux will use flux.has_soil = True
    Jv1 = - g.property('K')[v_base] * (psi_out[v_base] - psi_in[v_base])

    assert Jv1 < Jv3