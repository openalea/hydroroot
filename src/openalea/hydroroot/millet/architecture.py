# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 12:56:51 2015

@author: ndour

Revision 11/07/2025
- The module simulates the architecture of a millet root system
- It is quite similar of HydroRoot v1 architecture generator
- Q : collet seems to be an internode

"""
DEBUG = False

import random

import numpy as np

from openalea.mtg import *
from openalea.mtg import MTG, fat_mtg
from openalea.mtg import algo
from openalea.mtg.algo import orders, axis
from openalea.mtg.traversal import pre_order2_with_filter, post_order2

#from openalea.hydroroot import radius


# CPL: not used
# import numpy as np


def root_base(g):
    return next(g.component_roots_iter(g.root))


def add_branching_linear_axis(g, parent_id, nb_segments, **kwds):
    vertices = []
    vid = g.add_child(parent_id, edge_type='+', **kwds)
    vertices.append(vid)
    for i in range(nb_segments - 1):
        vid = g.add_child(vid, edge_type='<', **kwds)
        vertices.append(vid)
    return vertices


def millet_mtg(
    g=None,
    vid=None,
    seed=None,
    # seminal
    nb_vertices=350,
    branching_variability=0.2,
    branching_delay=0.2,
    LR_length_law=None,
    nude_tip_length=0,
    order_max=2,
    branching_stability=0.8,
    # crown
    nb_crown=2,
    crown_length= 50,
    segment_length=1e-4,
    **kwargs):
    """ Generate a Millet architecture.

    This parametric function generates a millet RSA.
    The set of parameters defines seminal and crown roots.

    The MTG will contain two type of roots :
        - Seminal
        - Crown

    Parameters
    ==========

        - g : MTG (None)
            The initial RSA
        - vid (int)
            vertex that will hold the generated sub-tree

    Returns
    =======

        - g : MTG

    Examples
    ========

        >>> g = millet_mtg()
    """

    nude_tip_length = (nb_vertices*15) / 75 # for a root 75cm-root we computed a 15cm nude tip length, so for a given root length we infer the corresponding length of apical zone
    # Review : this has to be set as input parameter
    #segment_length = 1e-4
    anchors = []  # list of Future branching points
    branching_positions = []  # list of distances to the apex of branching points
    collet_ids=[]

    if g is None:
        g = MTG()

    if vid is None :
        vid = g.add_component(g.root, label='collet', order = 0)
    collet_ids.append(vid)

    for i in range(nb_crown):
        vid= g.add_child(vid,label='collet',edge_type='<', order = 0)
        collet_ids.append(vid)


    def add_CR(nbre,length):
        label = 'Crown'

        for i in range(nbre):
            add_branching_linear_axis(g, collet_ids[i] , length, label=label)


    #if LR_length_law is None :
        #def LR_length_law(x):
            #return 10*np.random.geometric(0.77) # We use a geometric law to simulate LR length

    def delayed_markov(timer):
        """ markov chain with a delay between ramification """
        if (timer <= 0) :
            return (1,branching_delay)
        else :
            timer -= 1
            return 0,timer

    def create_randomized_delayed_axis(nid, n, anchors=anchors, **kwds):
        """ create an axis of length n using the delayed markov
            and randomized the id of the branching points in anchors
            around the theoretical branching positions

        :Parameters:
            - nid: root node for the axis
            - n : number of vertices for this axis
            - anchors: future ramification points on this axis
        """
        n = int(n)
        axis = []
        shuffled_axis = []
        branch, time = delayed_markov(0)


        for i in range(n-1):
            branch, time = delayed_markov(time)
            if (n-i) > (nude_tip_length): # we initialize the unbranched zone to zero
                axis.append((branch, n-i))  #n-i: distance to the tip
                shuffled_axis.append((branch,n-i))
            else:
                axis.append((0,0))
                shuffled_axis.append((0,0))

        for i in range(n-1):
            # shift (1-branching_stability) branching points
            # at max (1-branching_stability)*branching delay away from
            # theoretical branching position
            if (axis[i][0] == 1):   # read 'axis' only to avoid treating the same branching point after each shift
                if 1 :  #random.random() < branching_variability :
                    var = int(round(branching_variability*branching_delay))
                    shift = random.randint(-var,var)
                    # shift occurs only if the target is not branched already or outside the axis
                    if ((i+shift)>0) and ((i+shift)<n-1) and (shuffled_axis[i+shift][0]==0) : #the final position before the shift should not be branched
                        b, p = shuffled_axis[i] # we recover the mean branching point of branching point
                        shuffled_axis[i] = (0, p) #the branching point becomes an unbranched point
                        shuffled_axis[i+shift] = (1, p) # the new branching point after shift

        for ramif, position in shuffled_axis:
            # nid = g.add_child(nid, edge_type='<', **kwds)
            # nid_position_index = position
            # if ramif :
            #     anchors.append(nid)
            #     branching_positions.append(nid_position_index)
            order = nid.order
            nid = nid.add_child(order=order, edge_type='<',**kwds)
            nid.position_index = position
            if ramif :
                anchors.append(nid)
                #branching_positions.append(nid._vid)
        #print branching_positions

   

    if DEBUG:
        print('Begining MTG building')

    # Create the first axis
    create_randomized_delayed_axis(g.node(collet_ids[-1]), nb_vertices, label='Seminal')


    position_index = g.property('position_index')
    #Be sure that below a root length of 15cm, there is no branching
    if nb_vertices <= int(nude_tip_length):
        pass
    else:
        while anchors and (len(g) < nb_vertices*10):   # while they are branching point left
            # Review 11/07/25: we need a stop criteria
            nid = anchors.pop(0) 
            pos_index = nid.position_index # distance to the tip
            nid_order = nid.order
            if nid_order < order_max:  # check if maximal branching order was reached

                # if there is a length law, use it to compute lateral root length at this position
                if LR_length_law:
                        lateral_length = LR_length_law(pos_index)
                        
                else: 
                    # standard lateral root length - can't be longer than the bearing axis remaining branching length (remaining length - nude tip length)
                    # n = len(list(algo.descendants(g,nid._vid,RestrictedTo='SameAxis')))
                    # #n = random.randint(1, max(n-nude_tip_length,1))
                    # n = max(n-nude_tip_length,1)
                    lateral_length = int(50*np.random.geometric(0.77))


                # create axis of appropriate length
                if lateral_length > 0:
                    # branching_variability also apply to the length of LR
                    lateral_length = int(lateral_length)
                    # var = int(lateral_length*branching_variability)
                    # lateral_length = random.randint(max(1,lateral_length-var), lateral_length+var)
                    # Create the first  node of the branching point and the corresponding axis
                    cid = nid.add_child(order=nid_order+1, edge_type='+', label= g.label(nid._vid))
                    create_randomized_delayed_axis(cid,lateral_length, label=g.label(cid._vid))
                    #for i in range(lateral_length):
                        #cid= g.add_child(cid,edge_type='<')

    #print branching_positions
    #p_b = [position[i] for i in branching_positions]
    #t=zip(branching_positions,p_b)
    #print t
    if DEBUG:
        print('entering addition of crown root to the MTG')
    add_CR(nb_crown,crown_length)

    if DEBUG:
        print('ending crown roots addition to the MTG')

    if DEBUG:
        print('ending MTG building')

    g = fat_mtg(g)

    # Compute Order on each vertex
    _orders = orders(g, scale=g.max_scale())
    g.properties()['order'] = _orders

    return g


