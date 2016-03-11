import random

from openalea.mtg import MTG, fat_mtg
from openalea.mtg.traversal import post_order2
from openalea.mtg import algo


def mtg_builder(
    primary_length,
    primary_length_data,
    lateral_length_data,
    segment_length=1e-4,
    branching_variability=None,
    branching_delay=None,
    length_law=None,
    nude_tip_length=20,
    order_max=5,
    seed=None):
    """ Create a MTG from length laws.

    The first law is the length of the primary root.
    The second law is the length of the ramification.

    The primary root is of length 'total_length'.
    It is discretized following the law length_base.

    All the variables are expressed in meters.
    """
    length_base = primary_length_data
    length_lateral = lateral_length_data

    g = MTG()
    rid = vid = g.add_component(g.root)
    # Primary root

    # g.node(vid).tip_length = total_length
    prev_len = 0.
    ramifs = []

    for i, len_base in enumerate(length_base):
        while len_base-prev_len > 0:
            prev_len += segment_length
            vid = g.add_child(vid, edge_type='<', label='S', base_length=prev_len, length=segment_length, order=0)
        vid = g.add_child(vid, edge_type='<', label='S', base_length=len_base, length=segment_length, order=0)
        prev_len = len_base
        if length_lateral[i] > 0.:
            ramifs.append((length_lateral[i], vid))

    for length_ramif, rid in ramifs:
        len_tip = length_ramif
        cid = rid
        parent_base = g.node(cid).base_length
        prev_len = 0.
        while len_tip - prev_len >0.:
            prev_len += segment_length
            cid = g.add_child(cid, edge_type='+', label='S', base_length=prev_len+parent_base,
                            length=segment_length, order=1)

    ########################################################
    # Compute position_index: distance to tip in nb vertices
    position_index = g.properties()['position_index'] = {}
    for v in post_order2(g, rid):
        pi = position_index.setdefault(v,0)
        parent = g.parent(v)
        if parent is not None:
            if g.edge_type(v) == '<':
                position_index[parent]= pi+1


    ##################
    # Compute higher order of ramification
    anchors = [] # list of branching points where lateral roots will be "grafted"

    if not seed is None:
        random.seed(seed)

    def delayed_markov(timer):
        """ markov chain with a delay between ramification """
        if (timer <= 0) :
            return (1,branching_delay)
        else :
            timer -= 1
            return 0,timer

    def create_randomized_delayed_axis(nid, n, anchors=anchors):
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
            axis.append((branch, n-i))
            shuffled_axis.append((branch,n-i))

        for i in range(n-1):
            # shift (1-branching_stability) branching points
            # at max (1-branching_stability)*branching delay away from
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

    def update_randomized_delayed_axis(vid, anchors=anchors):
        """ create an axis of length n using the delayed markov
            and randomized the id of the branching points in anchors
            around the theoretical branching positions

        :Parameters:
            - nid: root node for the axis
            - n : number of vertices for this axis
            - anchors: future ramification points on this axis
        """
        n = len(list(algo.descendants(g,vid,RestrictedTo='SameAxis')))
        nid = g.node(vid)
        assert nid.order == 1

        axis = []
        shuffled_axis = []
        branch, time = delayed_markov(0)
        for i in range(n-1):
            branch, time = delayed_markov(time)
            axis.append((branch, n-i))
            shuffled_axis.append((branch,n-i))

        for i in range(n-1):
            # shift (1-branching_stability) branching points
            # at max (1-branching_stability)*branching delay away from
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
            if ramif :
                anchors.append(nid)

    # Update the ramification of order 1
    for v in g:
        if g.edge_type(v) == '+':
            update_randomized_delayed_axis(v)


    while anchors:   # while they are branching point left
        nid = anchors.pop(0)  # take next branching point
        position_index = nid.position_index # distance to the tip
        #print position_index
        if nid.order < order_max:  # check if maximal branching order was reached

            # if there is a length law, use it to compute lateral root length at this position
            if length_law:
                lateral_length = int(length_law(position_index))
            else : # standard lateral root length - can't be longer than the bearing axis remaining branching length (remaining length - nude tip length)
                n = len(list(algo.descendants(g,nid._vid,RestrictedTo='SameAxis')))
                #n = random.randint(1, max(n-nude_tip_length,1))
                n = max(n-nude_tip_length,1)
                lateral_length = n-1

            # create axis of appropriate length
            if lateral_length > 0:
                # branching_variability also apply to the length of LR
                var = int(lateral_length*branching_variability)
                lateral_length = random.randint(max(1,lateral_length-var), lateral_length+var)
                # Create the first  node of the branching point and the corresponding axis
                cid = nid.add_child(order=nid.order+1, edge_type='+')
                #print "pid length", nid, lateral_length
                create_randomized_delayed_axis(cid, lateral_length)

    g = fat_mtg(g)
    return g


