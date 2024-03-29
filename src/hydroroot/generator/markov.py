import random
import numpy as np

from binascii import hexlify as _hexlify
from os import urandom as _urandom

from openalea.mtg import *
from openalea.mtg import algo
from openalea.mtg import traversal
from hydroroot.law import length_law
#from random import choice

def my_seed():
    """generate my own seed function to capture the seed value."""
    return int(int(_hexlify(_urandom(2500)), 16) % 100000000)

def linear(n=5):
    """Create a MTG with just one linear axis without properties.

    :param n: (int) - number of vertices (Default value = 5)
    :returns: g (MTG)

    """
    g = MTG()
    root = g.root
    vid = g.add_component(root, order=0)
    for i in range(n - 1):
        vid = g.add_child(vid, edge_type='<', order=0)

    return g


def markov_binary_tree(g=None, vid=0, nb_vertices=1500,
                       branching_variability=0.1, branching_delay=20,
                       length_law=None,
                       nude_tip_length=200, order_max=5,
                       seed=None, censure_variability = False,  **kwargs):
    """
    Build a tree by modeling the stochastic process of lateral root branching by a first-order Markov chain on growing root axes.
    The building process is constrained by the arguments, see parameter descriptions.

    if `g` :code:`not None` then add the new tree to `g` at the vertex `vid`

    :param g: (MTG) (Default value = None)
    :param vid: (int) - the vertex id to which the new tree is added if not zero (Default value = 0)
    :param nb_vertices: (int) - the number of vertices in the main axis (Default value = 1500)
    :param branching_variability: (float) - variability on branching_delay and lateral lengths, :math:`\\in [0, 1]` (Default value = 0.1)
    :param branching_delay: (float) - number of vertex between branching (Default value = 20)
    :param length_law: function giving a length according to the branching position (Default value = None)
    :param nude_tip_length: (float) - distance from tip without branching in vertices number (Default value = 200)
    :param order_max: (int) - maximum root order (Default value = 5)
    :param seed: (int) - a seed for the random generators (Default value = None)
    :param censure_variability: (boolean) - deprecated because constrained by length_law- if True avoids to get a lateral longer than its parent (Default value = False)
    :returns:
        - g (MTG)

    .. seealso:: :func:`hydroroot.law.length_law`, :func:`hydroroot.law.histo_relative_law`, :func:`hydroroot.main.hydroroot_mtg`
    """
    # Modified FB 2020-03-10 : added flag in routine argument censure_variability see below
    if g is None:
        g = MTG()
    if vid == 0:
        vid = g.add_component(g.root, order=0)

    anchors = [] # list of branching points where lateral roots will be "grafted"

    if not seed is None:
        random.seed(seed)
        np.random.seed(seed)

    several_laws=True if isinstance(length_law, list) else False

    decrease = [ (1.-order/100.) for order in range(1, int(order_max)+1)]
    def markov():
        """simple random markov chain - unused now"""
        return 1 if random.random() < branching_variability else 0

    def delayed_markov(timer):
        """markov chain with a delay between ramification

        :param timer: 

        """
        if (timer <= 0) :
            return (1,branching_delay)
        else :
            timer -= 1
            return 0,timer

    def random_delayed_markov(timer):
        """random markov chain with a delay between
            possible ramification and uniform random variation
            of branching position around mean position

        :param timer: 

        """
        if (timer <= int(branching_variability*branching_delay)) :
            return (1,branching_delay) if (random.random() < (1-branching_variability)) else (0,timer)
        else :
            timer -= 1
            return 0,timer

    def create_axis(nid, n, anchors=anchors):
        """create a random axis of length n and record the id of the branching points in anchors

        :param nid: (int) - root node for the axis
        :param n: (int) - number of vertices for this axis
        :param anchors: (list) (Default value = anchors) - list of id of future ramification points on this axis

        """
        axis = [markov() for i in range(n)]
        for i in range(1,int(min(branching_delay,n))+1):
            axis[-i] = 0

        for ramif in axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            if ramif:
                anchors.append(nid)

    def create_delayed_axis(nid, n, anchors=anchors):
        """create an axis of length n using the delayed markov
            and record the id of the branching points in anchors

        :param nid: (int) - root node for the axis
        :param n: (int) - number of vertices for this axis
        :param anchors: (list) (Default value = anchors) - list of id of future ramification points on this axis

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

    def create_randomized_delayed_axis(nid, n, anchors=anchors):
        """create an axis of length n using the delayed markov
            and randomized the id of the branching points in anchors
            around the theoretical branching positions
        
        :param nid: (int) - root node for the axis
        :param n: (int) - number of vertices for this axis
        :param anchors: (list) (Default value = anchors) - list of id of future ramification points on this axis

        """
        n = int(n)
        axis = []
        shuffled_axis = []
        branch, time = delayed_markov(0)
        for i in range(n-1):
            branch, time = delayed_markov(time)
            axis.append((branch, n-i))
            shuffled_axis.append((branch,n-i))

        for i in range(n-1):
            # shift (1-branching_variability) branching points
            # at max (1-branching_variability)*branching delay away from
            # theoretical branching position
            if (axis[i][0] == 1):   # read 'axis' only to avoid treating the same branching point after each shift
                if 1 : #random.random() < branching_variability :
                    var = int(round(branching_variability*branching_delay))
                    shift = random.randint(-var,var)
                    # shift occurs only if the target is not branched already or outside the axis
                    if ((i+shift)>0) and ((i+shift)<n-1) and (shuffled_axis[i+shift][0]==0) :
                        b, p = shuffled_axis[i]
                        shuffled_axis[i] = (0, p)
                        shuffled_axis[i+shift] = (1, p)

        for ramif, position in shuffled_axis:
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<')
            nid.position_index = position
            if ramif :
                anchors.append(nid)

    # create_axis(g.node(vid), nb_vertices)  # deprecated

    #print 'entering MTG building'
    # Create the first axis
    create_randomized_delayed_axis(g.node(vid), nb_vertices)



    while anchors:   # while they are branching point left
        nid = anchors.pop(0)  # take next branching point
        position_index = nid.position_index # distance to the tip
        #print position_index
        current_order = nid.order
        if nid.order < order_max:  # check if maximal branching order was reached

            # if there is a length law, use it to compute lateral root length at this position
            if length_law:
                if several_laws:
                    current_law = length_law[0] if current_order == 0 else length_law[1]
                    lateral_length = int(current_law(position_index))
                else:
                    lateral_length = int(length_law(position_index))
            else : # standard lateral root length - can't be longer than the bearing axis remaining branching length (remaining length - nude tip length)
                n = len(list(algo.descendants(g,nid._vid,RestrictedTo='SameAxis')))
                #n = random.randint(1, max(n-nude_tip_length,1))
                n = max(n-nude_tip_length,1)
                lateral_length = n-1

            real_lateral_length = lateral_length
            # create axis of appropriate length
            if lateral_length > 0:
                # branching_variability also apply to the length of LR
                var = int(lateral_length*branching_variability)
                lateral_length = random.randint(max(1,lateral_length-var), lateral_length+var)

                # Censure variability
                # Modified FB 2020-03-10 : by default lateral lengths not constrained to nb_vertices*decrease[current_order]
                #       so it is an "arbitrary" constrain on length
                # the lateral lengths are constrained to the length law in def histo_relative_law
                if censure_variability:
                    if lateral_length > nb_vertices*decrease[current_order]:
                        if real_lateral_length <= nb_vertices*decrease[current_order]:
                            lateral_length = real_lateral_length
                        else:
                            print(("WARNING: lateral length is too large ", lateral_length))
                            lateral_length = nb_vertices*decrease[current_order]



                # Create the first  node of the branching point and the corresponding axis
                cid = nid.add_child(order=nid.order+1, edge_type='+')
                #print "pid length", nid, lateral_length
                #print "Var, Lateral length: ", var, lateral_length

                create_randomized_delayed_axis(cid, lateral_length)

    fat_mtg(g)
    #print 'exiting MTG building'
    return g


def shuffle_axis(g=None, shuffle=False):
    """For each subtree of a MTG, change its root node to another node of the same axis.

    :param g:  (Default value = None)
    :param shuffle:  (Default value = False)
    :returns: - g

    """
    max_scale = g.max_scale()
    if shuffle:
        #print 'entering axis shuffling'
        v_base = next(g.component_roots_at_scale_iter(g.root, scale=max_scale))
        shuffling = {}
        shuffled = []

        ramifs = (v for v in g.vertices_iter(scale=max_scale) if g.edge_type(v) == '+')

        for v in ramifs:
            axis = g.Axis(g.parent(v))    # list of all node in the same axis as v
            shuffling[v] = random.choice(axis)  # record new position for each subtree

        for v in traversal.post_order2(g,v_base):
            if v in list(shuffling.keys()) and v not in shuffled:
                sub = g.sub_tree(v, copy=True)        # get subtree
                g.remove_tree(v)                      # prune it from old position
                newbranch = g.add_child_tree(shuffling[v], sub)  # insert subtree at previously determined position
                shuffled.append(newbranch[0])        # keep track of shifted tree id
        #print 'exiting axis shuffling'
    return g


def generate_g(seed = None, length_data = None, branching_variability = 0.25,
               delta = 2e-3, nude_length = 2e-3, primary_length = 0.13, segment_length = 1e-4, order_max = 4):
    """generate a MTG according to the input parameters using length_data (mendatory) for the branching

    Preprocessing variables before calling :func:`markov_binary_tree`:

        - nb_vertices: primary root length in vertex number
        - branching_delay: distance between branching in vertex number
        - nb_nude_vertices: distance from tip without branching in vertex number
        - length_law: the function giving the lateral length according to its position on the parent branch
    
    :param seed: (int) - the seed for the random generator in the markof chain (Default value = None)
    :param length_data: (Dataframe) - pandas dataframe columns names 'LR_length_mm' and 'relative_distance_to_tip' sorted by 'relative_distance_to_tip' (Default value = None)
    :param branching_variability: (float) - variability on branching_delay and lateral lengths, :math:`\\in [0, 1]` (Default value = 0.25)
    :param delta: (float) - reference distance (m) between successive branching axis (Default value = 2e-3)
    :param nude_length: (float) - length (m) at root tip with no ramification (Default value = 2e-3)
    :param primary_length: (float) - primary root length (m) (Default value = 0.13)
    :param segment_length: (float) - length of the vertices (m) (Default value = 1e-4)
    :param order_max: (int) - maximum lateral roots order (Default value = 4)
    
    :returns:
        - g: MTG with the following properties set: edge_type, label, position

    .. seealso:: :func:`hydroroot.law.length_law`, :func:`hydroroot.law.histo_relative_law`, :func:`hydroroot.main.root_builder`
    """

    # nude length and branching delay in terms of number of vertices
    nb_nude_vertices = int(nude_length / segment_length)
    branching_delay = int(delta / segment_length)

    nb_vertices = int(primary_length / segment_length)

    # F. Bauget 2022-08-12: added if-else to be able to use the function without length data, useful for usage demo
    if length_data:
        length_max_secondary = length_data[0].LR_length_mm.max() * 1e-3  # in m

        law_order1 = length_law(length_data[0], scale_x = primary_length / 100., scale = segment_length)
        law_order2 = length_law(length_data[1], scale_x = length_max_secondary / 100., scale = segment_length)
        _length_law = [law_order1, law_order2]
    else:
        _length_law = None

    g = markov_binary_tree(
        nb_vertices = nb_vertices,
        branching_variability = branching_variability,
        branching_delay = branching_delay,
        length_law = _length_law,
        nude_tip_length = nb_nude_vertices,
        order_max = order_max,
        seed = seed)
    return g
