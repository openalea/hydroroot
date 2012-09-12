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
                       branching_chance=0.1, branching_delay=50, 
                       length_law=None,
                       nude_tip_length=200,  order_max=5, 
                       seed=None,  **kwargs ):
    """
    Parameters
    ----------
        - g : MTG
        - vid : id of the root of the MTG where the generating tree will be added
        - nb_vertices : number of element of the main axis
        - branching_chance : probability of ramification at each point
        - branching_delay : reference distance between successive branching axis
        - length_law : spline given the length of lateral ramification
        - nude_tip_length : length at root tip with no ramification
        - LR_length_law : distribution of LR length along axis length
        - seed : Seed for random number generator (default=None).
    """
    if g is None:
        g = MTG()
    if vid == 0:
        vid = g.add_component(g.root, order=0)

    anchors = [] # list of branching points where lateral roots will be "grafted"

    if not seed is None:
        random.seed(seed)

    def markov():   
        """ simple random markov chain - unused now """
        return 1 if random.random() < branching_chance else 0

    def delayed_markov(timer):    # random markov chain with threshold and a delay between possible ramification
        if (timer == 0) :
            return (1,branching_delay) if (random.random() < branching_chance) else (0,0)
        else :
            timer -= 1
            return 0,timer

    def create_axis(nid, n, anchors=anchors):    
        """ create a random axis of length n and record the id of the branching points in anchors """
        axis = [markov() for i in range(n)]
        for i in range(1,min(branching_delay,n)+1):
            axis[-i] = 0

        for ramif in axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            if ramif:
                anchors.append(nid)

    def create_delayed_axis(nid, n, anchors=anchors):   
        """ create an axis of length n using the delayed markov 
            and record the id of the branching points in anchors

        :Parameters:
            - nid: root node for the axis
            - n : number of vertices for this axis
            - anchors: future ramification points on this axis
        """
        axis = []
        branch, time = delayed_markov(0)
        for i in range(n-1):
            branch, time = delayed_markov(time)
            if (n-i) > (nude_tip_length): # check axis length compared to minimal branching length
                axis.append((branch, n-i))
            else : # leave end of axis empty of branching
                axis.append((0,0))
        for ramif, position in axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            nid.position_index = position
            if ramif:
                anchors.append(nid)

    # create_axis(g.node(vid), nb_vertices)  # deprecated

    # Create the first axis
    create_delayed_axis(g.node(vid), nb_vertices)

    while anchors:   # while they are branching point left
        nid = anchors.pop(0)  # take next branching point
        position_index = nid.position_index # distance to the tip
        if nid.order < order_max:  # check if maximal branching order was reached
            # Create the first node of the branching point
            cid = nid.add_child(order=nid.order+1, edge_type='+')
            # Compute length of root downward of the branching point
            n = len(list(algo.descendants(g,nid._vid,RestrictedTo='SameAxis')))
            # New lateral root can't be longer than the bearing axis remaining branching length (remaining length - nude tip length)
            n = random.randint(1, max(n-nude_tip_length,1))
            # Create the new axis
            #create_axis(cid, n-1) # deprecated
            lateral_length = n-1
            if length_law:
                lateral_length = int(length_law(position_index))

            if lateral_length:
                create_delayed_axis(cid, lateral_length)

    fat_mtg(g)
    return g




