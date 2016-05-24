from hydroroot.main import hydroroot as hydro


def data():
    length = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
    axial = ([0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18],
        [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4, 292.6e-4, 262.5e-4, 511.1e-4])
    radial = ([0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16],
        [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300])

    return length, axial, radial

def closed(delta, eps=1e-14):
    assert(abs(delta) < eps)

def check_flux(g, Jv_global):

    # Jv_global is the sum of radial fluxes
    j = g.property('j')
    delta = sum(j.itervalues()) - Jv_global

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

def test_flux1():
    length, axial, radial = data()
    g, surface, volume, Keq, Jv_global = hydro(primary_length=0.09,
                                               order_decrease_factor=0.7,
                                               length_data=length,
                                               axial_conductivity_data=axial,
                                               radial_conductivity_data=radial,
                                               seed=2)

    assert len(g) == 8584
    assert(3.7e-4 < surface < 3.8e-4)
    assert(1.3e-8 < volume < 1.4e-8)
    assert(0.1029 < Keq < 0.103)
    assert(0.03 < Jv_global < 0.031)

    check_flux(g, Jv_global)

    rid = 1
    lid = 899
    assert(g.is_leaf(lid))
    assert(g.parent(rid) is None)
    mid = 450

    psi_out = g.property('psi_out')

    closed(psi_out[rid]-0.1)
    closed(psi_out[lid]-0.15145569239681866)
    closed(psi_out[mid]-0.13587579341647765)

    # Test if Keq increase on the Trunk
    #check_Keq_monotony(g)
    return g, surface, volume, Keq, Jv_global


