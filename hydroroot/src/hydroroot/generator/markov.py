import random

from openalea.mtg import *
from openalea.mtg import algo

def linear(n=5):
    """ 
    Create a MTG with just one linear axis without properties.
    """
    g = MTG()
    root = g.root
    vid = g.add_component(root, order=0)
    for i in range(n-1):
        vid = g.add_child(vid, edge_type='<', order=0)

    return g

def markov_binary_tree(g=None, vid=0, nb_vertices=300,
                       branching_chance=0.1, branching_delay=50, order_max=5, seed=None,  **kwargs ):
    """
    Parameters
    ----------
        - g : MTG
        - vid : id of the root of the MTG where the generating tree will be added
        - nb_vertices : number of element of the main axis
        - branching_chance : probability of ramification at each point
        - branching_delay: minimal distance between the tip and the first branching axis 
        - seed : Seed for random number generator (default=None). 
    """
    if g is None:
        g = MTG()
    if vid == 0:
        vid = g.add_component(g.root, order=0)

    anchors = [] # branching points

    if not seed is None:
        random.seed(seed)

    # First axis
    def markov():
        return 1 if random.random() < branching_chance else 0

    def delayed_markov(timer):
        if (timer == 0) :
            return (1,branching_delay) if (random.random() < branching_chance) else (0,0)
        else :
            timer -= 1
            return 0,timer

    def create_axis(nid, n, anchors=anchors):
        axis = [markov() for i in range(n)]
        for i in range(1,min(branching_delay,n)+1):
            axis[-i] = 0

        for ramif in axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            if ramif:
                anchors.append(nid)

    def create_delayed_axis(nid, n, anchors=anchors):
        axis = []
        branch, time = delayed_markov(0)
        for i in range(n-1):
            branch, time = delayed_markov(time)
            if (n-i) > (100):
                axis.append(branch)
            else : # leave end of axis empty of branching
                axis.append(0)
        for ramif in axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            if ramif:
                anchors.append(nid)

    #create_axis(g.node(vid), nb_vertices)
    create_delayed_axis(g.node(vid), nb_vertices)

    while anchors:
        nid = anchors.pop(0)
        if nid.order < order_max:
            # Create a new axis
            cid = nid.add_child(order=nid.order+1, edge_type='+')
            n = len(list(algo.descendants(g,nid._vid,RestrictedTo='SameAxis')))
            n = random.randint(1, n)

            #create_axis(cid, n-1)
            create_delayed_axis(cid, n-1)

    fat_mtg(g)
    return g


def test_linear(n=30, psi_e=0.3, psi_base=0.1, Jv=15, k0=0.5, p_cst = 50.):
    """ Test flux and water potential computation on a linear root. """
    # topology
    #n=40

    g = linear(n)
    g = compute_flux(g,n=n, psi_e=psi_e, psi_base=psi_base, Jv=Jv, k0=k0, p_cst = p_cst)
    scene = plot(g, prop_cmap='J_out', has_radius=True)
    return g, scene


def test_tree(n=30,psi_e=0.3, psi_base=0.1, Jv=15, k0=0.5, p_cst = 50.):
    """ Test flux and water potential computation on a linear root. """
    # topology
    #n=40

    g = markov_binary_tree(nb_vertices=n)
    g = compute_flux(g,n=n, psi_e=psi_e, psi_base=psi_base, Jv=Jv, k0=k0, p_cst = p_cst)
    scene = plot(g, prop_cmap='J_out', has_radius=True)
    return g, scene


