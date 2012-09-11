# License

"""

.. todo:: Compute the number of xylem pipes in each root segment
and the radius of each xylem pipe.

"""
from openalea.mtg import *
from openalea.mtg import algo


def cont_radius(g, r_base, r_tip):
    """ Compute the radius of each segment of a root system.

    Set radius for elements of a mtg with an increase rate computed from 
    given base and tip radius in a continuous way.
    """
    r_base, r_tip = float(r_base), float(r_tip)

    assert (r_base>r_tip),"Base radius should be greater than tip radius"

    base = g.component_roots(g.root).next()
    base = g.node(base)
    base.radius = r_base

    _tips = dict((vid, g.order(vid)) for vid in g.vertices(scale = g.max_scale()) if g.is_leaf(vid))
    tips = {}
    for tip,order in _tips.iteritems():
        tips.setdefault(order, []).append(tip)

    max_order = max(tips)
    for order in range(max_order+1):
        for tip in tips[order]:
            l = [g.node(vid) for vid in algo.axis(g,tip)]
            n = len(l)
            if n==1:  # only one segment in the lateral root - very young lateral
                l[0].radius = r_tip
            else :    # more than one segment in the lateral root
                parent = l[0].parent()
                r0 = parent.radius if parent else r_base
                dr = (r0-r_tip)/(n-1)
                for node in l:
                    if node.radius==r0:  # keep the value of the base/junction radius for the first segment  of the lateral
                        continue
                    node.radius = node.parent().radius - dr  # decrease the radius from base to tip

    return g

def discont_radius(g, r_base, r_tip):
    """ Compute the radius of each segment of a root system.

    Set radius for elements of a mtg with an increase rate computed
    from the length of the longest axis and its base and tip radius.
    
    Radius can be discontinuous e.g. for a young/small lateral on an old root,
    the yound root radius is very small initially compared to the old one.
    """
    r_base, r_tip = float(r_base), float(r_tip)

    assert (r_base>r_tip),"Base radius should be greater than tip radius"

    base = g.component_roots(g.root).next()
    base = g.node(base)
    base.radius = r_base

    _tips = dict((vid, g.order(vid)) for vid in g.vertices(scale = g.max_scale()) if not algo.sons(g,vid,EdgeType='<'))
    tips = {}
    for tip,order in _tips.iteritems():
        tips.setdefault(order, []).append(tip)

    max_order = max(tips)

    max_len = 0
    for order in range(max_order+1):   #find the longest axis length among all axis
        for tip in tips[order]:
            max_len = max(max_len, len(list(algo.axis(g,tip))))
    assert (max_len>1), "MTG too short for analysis"
    growth_rate = (r_base-r_tip)/(max_len-1)    #define growth rate according to radius extremities of the longest axis

    # radius are computed from tips to bases according to growth rate extrapolated from absolute longest axis of the MTG
    for order in range(max_order+1):
        for tip in tips[order]:
            node = g.node(tip)
            node.radius = r_tip
            while node and node.parent() and node.edge_type() != '+':
                node.parent().radius = node.radius + growth_rate
                node = node.parent()

    return g

def compute_length(g, length = 1.e-4):
    length = float(length)
    for vid in g.vertices(scale=g.max_scale()):
        g.node(vid).length = length
    return g

def compute_relative_position(g):
    scale = g.max_scale()
    position = {}
    position_measure = {}
    axis_length = {}
    root_id = g.component_roots_at_scale(g.root, scale=scale).next()
    for vid in traversal.post_order2(g, root_id):
        sons = algo.sons(g,vid,EdgeType='<')
        position[vid] = position[sons[0]]+1 if sons else 0
        position_measure[vid] = position[vid]*g.node(vid).length
        if g.edge_type(vid) == '+' or g.parent(vid) is None:
            axis_length[vid] = position[vid]

    relative_position = {}
    for axis_id, _length in axis_length.iteritems():
        _length = float(max(1, _length))
        for v in algo.local_axis(g,axis_id):
            relative_position[v] = position[v] / _length

    g.properties()['position'] = position_measure
    g.properties()['relative_position'] = relative_position

    return g

