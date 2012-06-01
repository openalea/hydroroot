import random

from hydroroot.flux import *
from openalea.mtg import *
from openalea.mtg import algo
from openalea.mtg import turtle as turt

def tap_root(n=5):
    """ 
    Create a MTG with just one linear axis without properties.
    """
    g = MTG()
    root = g.root
    vid = g.add_component(root)
    for i in range(n-1):
        vid = g.add_child(vid, edge_type='<')

    return g


#BUG : the produced MTG do not have node 1 ?
def random_branched_root(n_elements=20, branching_chance=0.25, max_branch_height=5, branching_delay = 3):
    """ Create a random MTG of max n elements with a given chance to branch 
    (no multiple branching per element allowed, 
    no branching on an axis until the branching_delay is past) 
    and random branch length (up to max_branch_height).
    """
    g = MTG()
    root = g.root
    vid = g.add_component(root)

    n = int(n_elements)
    order = 1
    branch_height = {}    # keep track of branch length still left to build
    branching_point = {}  # keep track of branching points id
    branch_elements = n   # keep track of the number of elements still available to build the MTG
    delay = branching_delay   # keep track of the delay between successive branching events

    # budget the length of the main axis
    branch_height[1] = random.randint(max(1,delay),min(max_branch_height,branch_elements))
    branch_elements -= branch_height[1]

    # build the MTG
    for i in range(n-1):
        print order,branch_height[order]
        print "left", branch_elements

        # branching only if there are elements left to build new branches & 
        #                   there is no branch already at the current node & 
        #                   the branching delay is passed & 
        #                   we are not too close from the tip (more than branching_delay element left in the supporting axis)
        if ((random.random()<branching_chance) and 
            (branch_elements) and 
            (g.node(vid).nb_children() == 0) and 
            (delay == 0) and 
            (branch_height[order]>branching_delay)):             

            order += 1   # go up one scale at each branching point
            delay = branching_delay   # block future branching for a given time
            branching_point[order] = vid   # keep track of the branching point id at the current scale to come back there at the end of the branch
            branch_height[order] = random.randint(1,min(max_branch_height,branch_elements)) # new branch can use up to max_branch_height elements or what is left over from the pool of elements
            branch_elements -= branch_height[order]  # keep track of the number of elements left to build the rest of the MTG
            vid = g.add_child(vid, edge_type='+', scale_id = order)  # add the branch element
            branch_height[order] -= 1 # decrease the length of the branch left to build at this scale
            if not branch_height[order] : # last branch element was added, end the branch
                vid = branching_point[order]  # next element will be added as a child to the last branching point
                order -= 1
        elif (branch_height[order]):   # no branching, add child element linearly if there are still some left to add in the current axis
            vid = g.add_child(vid, edge_type='<', scale_id = order)
            branch_height[order]-= 1  # decrease length of axis left to build
            if (order>1) and not (branch_height[order]) : # if it was the last element of a branch, get down one scale
                vid = branching_point[order] # next element will be added as a child to the last recorded branching point
                order -= 1
            if delay :  # decrease delay before the next branching if it is not passed already
                delay -= 1

    return g



#BUG : the produced MTG do not have node 1 ?
def grown_branched_root(n_elements=20, branching_chance=0.25, branching_delay = 3):
    """ Create a MTG of n elements by growing it with a given chance to branch (no multiple branching per element allowed, no branching on an axis until the branching_delay is past).

    """
    g = MTG()
    root = g.root
    vid = g.add_component(root)

    n = int(n_elements)
    branch_count = 1                     # main root is labelled 1
    branches_tips = {}                   # keep track of branches tips
    branches_tips[branch_count] = vid    # initial tip of the main root
    branch_elements = n                  # keep track of the number of elements still available to build the MTG
    delay = {}                             # keep track of branching clock of the different axis
    delay[branch_count] = branching_delay   # initial clock for the first axis

    # build the MTG
    while (branch_elements>0) :

        for branch in branches_tips.keys() :  # for each axis, either grow or branch
            reserved = branching_delay*branch_count    # keep enough elements to end each axis by a non branched tip
            if (random.random()<branching_chance) and (branch_elements>reserved) and (g.node(branches_tips[branch]).nb_children() == 0) and (delay[branch] == 0) :   # branching only if they are elements left to build new branches & there is no branch already at the current node & the branching delay is passed & we are not too close from the tip (branching_delay elements left for each axis)
                delay[branch] = branching_delay  # reset clock for the bearing branch
                branch_count += 1   #  increase branch count
                delay[branch_count] = branching_delay  # start the clock for the new branch
                vid = g.add_child(branches_tips[branch], edge_type='+', branch_id = branch_count)  # add first element of new branch
                branches_tips[branch_count] = vid    # record the id of new branch tip
                branch_elements -= 1   # decrease the number of available elements
            elif (branch_elements>0) :   # grow
                vid = g.add_child(branches_tips[branch], edge_type='<', branch_id = branch)   # add element to existing branch
                branches_tips[branch] = vid  # record this element as the new tip
                branch_elements -= 1   # decrease the number of available elements
                if (delay[branch]>0) :  # decrease branching clock for current branch
                    delay[branch] -= 1
    return g

def random_binary_tree(g=None, vid=0, nb_vertices=20, 
                       branching_chance=0.25, branching_delay=3, order_max=5, **kwargs ):
    """
    Parameters
    ----------
        - g : MTG
        - vid : id of the root of the MTG where the generating tree will be added
        - nb_vertices : number of element of the main axis
        - branching_chance : probability of ramification at each point
        - branching_delay: minimal distance between the tip and the first branching axis 
    """
    if g is None:
        g = MTG()
    if vid == 0:
        vid = g.add_component(g.root, order=0)

    anchors = [] # branching points
    
    # First axis
    def markov():
        return 1 if random.random() < branching_chance else 0

    def create_axis(nid, n, anchors=anchors):
        axis = [markov() for i in range(n)]
        for i in range(1,min(branching_delay,n)+1):
            axis[-i] = 0

        for ramif in axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            if ramif: 
                anchors.append(nid)
    
    create_axis(g.node(vid), nb_vertices)

    while anchors:
        nid = anchors.pop(0)
        if nid.order < order_max:
            # Create a new axis
            cid = nid.add_child(order=nid.order+1, edge_type='+')
            n = len(list(algo.descendants(g,nid._vid,RestrictedTo='SameAxis')))
            n = random.randint(1, n)

            create_axis(cid, n-1)
    fat_mtg(g)
    return g

def cont_radius(g, r_base, r_tip):
    """ Set radius for elements of a mtg with an increase rate computed from given base and tip radius in a continuous way.
    """

    assert (r_base>r_tip),"Base radius should be greater than tip radius"

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


def discont_radius(g, r_base, r_tip):
    """ Set radius for elements of a mtg with an increase rate computed from the length of the longest axis and its base and tip radius.
        Radius can be discontinuous e.g. for a young/small lateral on an old root, the yound root radius is very small initially compared to the old one.
    """

    assert (r_base>r_tip),"Base radius should be greater than tip radius"

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
    
    max_len = 0
    for order in range(max_order+1):   #find the longest axis length among all axis
        for tip in tips[order]:
            max_len = max(max_len, len(list(algo.axis(g,tip))))
    assert (max_len>1), "MTG too short for analysis"
    growth_rate = (r_base-r_tip)/(max_len-1)    #define growth rate according to radius extremities of the longest axis

    # radius are computed from tips to bases according to growth rate extrapolated from absolute longest axis of the MTG
    for order in range(max_order+1):
        for tip in tips[order]:
            l = [g.node(vid) for vid in algo.axis(g,tip)]
            for i in range(len(l)):
                node = l.pop()   # this assume that the l list is ordered, with root tip at the last place
                if g.is_leaf(node.index()):   # root tip have small radius
                    node.radius = r_tip
                else :  # root segment get bigger as we get away from the tip
                    suc_index = algo.successor(g,node.index())
                    node.radius = g.node(suc_index).radius + growth_rate


def compute_k(g, k0 = 0.1):
    """ Set radial conductances in a MTG at a given value. """
    k = dict.fromkeys(g.vertices(scale=g.max_scale()), k0)
    return k


def compute_K(g, cte=10):
    """ Set axial conductances in a MTG according to Poiseuille law. """
    def poiseuille(r):
        return cte*(r**4)
    radius = g.property('radius')
    K = dict((vid, poiseuille(radius[vid])) for vid in g.vertices(scale=g.max_scale()))
    return K


def test_linear(n=5, psi_e=0.3, psi_base=0.1, Jv=15, k0=0.5, p_cst = 50.):
    """ Test flux and water potential computation on a linear root. """
    # topology
    #n=40

    k0 = float(k0) / n
    p_cst = float(p_cst) / n
    g = tap_root(n)
    discont_radius(g, r_base=10, r_tip=1)
    k = compute_k(g, k0)
    K = compute_K(g,p_cst)

    assert all(v>0 for v in K.values()),K

    g = flux(g, k, K, Jv, psi_e, psi_base)

    J_out = g.property('J_out')
    assert all(v>0 for v in J_out.values()), J_out.values()
    return g


def plot(length=5., r_base=1., r_tip=0.25):
    """
    Exemple:

        >>> from openalea.plantgl.all import *
        >>> s = plot()
        >>> shapes = dict( (x.getId(), x.geometry) for x in s)
        >>> Viewer.display(s)
    """
    g = random_binary_tree()
    print 'Nb vertices ', len(g)
    # compute length
    for v in g:
        n = g.node(v)
        n.length = length 
    
    discont_radius(g,r_base=r_base, r_tip=r_tip)

    angles = [90,45]+[30]*5
    def root_visitor(g, v, turtle):
        n = g.node(v)
        radius = n.radius
        order = n.order
        length = n.length

        if g.edge_type(v) == '+':
            angle = angles[order]
            turtle.down(angle)


        turtle.setId(v)
        turtle.setWidth(radius)
        for c in n.children():
            if c.edge_type() == '+':
                turtle.rollL(180)
        turtle.setColor(order+1)
        turtle.F(length)

    turtle = turt.PglTurtle()
    turtle.down(180)
    scene = turt.TurtleFrame(g, visitor=root_visitor, turtle=turtle, gc=False)
    return scene
