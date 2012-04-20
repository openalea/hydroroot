from hydroroot.flux import *
from openalea.mtg import *
from openalea.mtg import algo

def tap_root(n=5):
    """ Create a MTG with just one linear axis without properties."""
    g = MTG()
    root = g.root
    vid = g.add_component(root)
    for i in range(n-1):
        vid = g.add_child(vid, edge_type='<')

    return g

def radius(g, r_base, r_tip):
    base = g.component_roots(g.root).next()
    base = g.node(base)
    base.radius = r_base

    _tips = dict((vid, g.order(vid)) for vid in g.vertices(scale = g.max_scale()) if g.is_leaf(vid))
    tips = {}
    for tip,order in _tips.iteritems():
        tips.setdefault(order, []).append(tip)
#X         Same code than setdefaults
#X         if order in tips:
#X             tips[order].append(tip)
#X         else:
#X             tips[order] = [tip]

    max_order = max(tips)
    for order in range(max_order+1):
        for tip in tips[order]:
            l = [g.node(vid) for vid in algo.axis(g,tip)]
            n = len(l)
            parent = l[0].parent()
            r0 = parent.radius if parent else r_base
            dr = (r_tip-r0)/n
            for node in l:
                if node.radius:
                    continue
                node.radius = node.parent().radius + dr



def compute_k(g, k0 = 0.1):
    k = dict.fromkeys(g.vertices(scale=g.max_scale()), k0)
    return k

def compute_K(g, cte=10):

    def poiseuille(r):
        return cte*(r**4)

    radius = g.property('radius')
    K = dict((vid, poiseuille(radius[vid])) for vid in g.vertices(scale=g.max_scale()))
    return K

def test_linear(n=5, psi_e=100, psi_base=1, Jv=100, k0=0.5, p_cst = 50.):
    # topology
    #n=40
    
    k0 = float(k0) / n
    p_cst = float(p_cst) / n
    g = tap_root(n)
    radius(g, r_base=1, r_tip=1)
    k = compute_k(g, k0)
    K = compute_K(g,p_cst)

    assert all(v>0 for v in K.values()),K

    g = flux(g, k, K, Jv, psi_e, psi_base)

    J_out = g.property('J_out')
    assert all(v>0 for v in J_out.values()), J_out.values()
    return g
