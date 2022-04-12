from hydroroot.main import hydroroot as hydro
from hydroroot import flux


def data():
    length = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
    axial = ([0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18],
        [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4, 292.6e-4, 262.5e-4, 511.1e-4])
    radial = ([0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16],
        [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300])

    return length, axial, radial

def closed(delta, eps=1e-14, txt=''):
    if not txt:
        assert abs(delta) < eps
    else:
        assert abs(delta) < eps, txt


def check_flux(g, Jv_global):

    # Jv_global is the sum of radial fluxes
    j = g.property('j')
    delta = sum(j.values()) - Jv_global

    closed(delta)

def check_Keq_monotony(g):
    """
    Error the Keq is increasing while it would normally decrease.
    """
    Keq = g.property('Keq')
    trunk = g.Trunk(1)
    n = len(trunk)

    keq_prev = Keq[trunk[0]]
    for i in range(1, n):

        assert -1e-16 <= keq_prev - Keq[trunk[i]] , 'Failure : vertex %d, value: %f' %(trunk[i], Keq[trunk[i]] < keq_prev)
        keq_prev = Keq[trunk[i]]

def check_radial_flow_conservation(g):
    """j = ki *(psi_e -psi_in)

    """
    psi_in = g.property('psi_in')
    psi_e = 0.4
    k = g.property('k')
    j = g.property('j')

    for vid in j:
        ji = j[vid]
        ki = k[vid]
        psi_i = psi_in[vid]

        delta = ji - ki*(psi_e-psi_i)

        closed(delta, txt = "For vertex %d, radial flow is not propoerly calculated."%vid)

def check_vertex_flow_conservation(g):
    """j = ki *(psi_e -psi_in)

    """
    J_out = g.property('J_out')
    j = g.property('j')

    for vid in j:

        J = J_out[vid]
        J_c = 0.
        for cid in g.children(vid):
            J_c += J_out[cid]

        ji = j[vid]
        delta = J - (ji + J_c)

        closed(delta, txt = "For vertex %d, Millmans law is not satisfied."%vid)



def test_flux1():
    length, axial, radial = data()
    g, surface, volume, Keq, Jv_global = hydro(primary_length=0.09,
                                               order_decrease_factor=0.7,
                                               length_data=length,
                                               axial_conductivity_data=axial,
                                               radial_conductivity_data=radial,
                                               seed=2)
    # F. Bauget 2022-04-11: the markov_binary_tree has been changed so the values were not correct anymore
    assert len(g) == 9896
    assert(4.1e-4 < surface < 4.2e-4)
    assert(1.4e-8 < volume < 1.5e-8)
    assert(0.1155 < Keq < 0.1156)
    assert(0.034 < Jv_global < 0.035)

    check_flux(g, Jv_global)

    rid = 1
    lid = 899
    assert(g.is_leaf(lid)) #no leaf
    assert(g.parent(rid) is None)
    mid = 450

    psi_out = g.property('psi_out')

    closed(psi_out[rid]-0.1)
    closed(psi_out[lid]-0.15120969955387603)
    closed(psi_out[mid]-0.1357112706017493)

    # Test if Keq increase on the Trunk
    #check_Keq_monotony(g)

    check_radial_flow_conservation(g)
    check_vertex_flow_conservation(g)

    return g, surface, volume, Keq, Jv_global


def test_cut():
    length, axial, radial = data()
    g, surface, volume, Keq, Jv_global = hydro(primary_length=0.09,
                                               order_decrease_factor=0.7,
                                               length_data=length,
                                               axial_conductivity_data=axial,
                                               radial_conductivity_data=radial,
                                               seed=2)

    g_cut = flux.cut(g, 0.04, threshold = 1e-4) # g_cut = flux.cut(g, 0.04)
    check_length(g_cut, 0.04, segment_length=1e-4)

def test_cut_and_flow():
    length, axial, radial = data()
    g, surface, volume, Keq, Jv_global = hydro(primary_length=0.09,
                                               order_decrease_factor=0.7,
                                               length_data=length,
                                               axial_conductivity_data=axial,
                                               radial_conductivity_data=radial,
                                               seed=2)

    # g_cut = flux.cut(g, 0.04, threshold = 1e-4) #g_cut = flux.cut(g, 0.04)
    g_cut = flux.cut_and_set_conductance(g, 0.04, threshold = 1e-4)
    check_length(g_cut, 0.04, segment_length=1e-4)

    g_cut = flux.flux(g_cut, cut_and_flow=True, invert_model=True)
    v_base = 1
    psi_e = 0.4
    psi_base = 0.101325
    Keqs = g_cut.property('Keq')
    Jv_cut = Keqs[v_base] * (psi_e - psi_base)

    print(Jv_cut)
    assert Jv_cut > Jv_global

def check_length(g, length, segment_length):
    axis = g.Axis(1)
    n = len(axis)
    _length = n*segment_length

    assert abs(_length-length) <= segment_length

