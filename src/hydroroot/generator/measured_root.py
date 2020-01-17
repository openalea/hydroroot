import random
import numpy as np

from openalea.mtg import MTG, fat_mtg
from openalea.mtg.traversal import post_order2
from openalea.mtg import algo

SUPERIOR_ORDER = True

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

    # Create the primary axis
    for i, len_base in enumerate(length_base):
        while len_base-prev_len > 0:
            prev_len += segment_length
            vid = g.add_child(vid, edge_type='<', label='S', base_length=prev_len, length=segment_length, order=0)
        vid = g.add_child(vid, edge_type='<', label='S', base_length=len_base, length=segment_length, order=0)
        prev_len = len_base
        if length_lateral[i] > 0.:
            ramifs.append((length_lateral[i], vid))

    # Create lateral axis
    for length_ramif, cid in ramifs:
        len_tip = length_ramif
        # print 'length ', len_tip
        _root_id = cid
        parent_base = g.node(cid).base_length
        prev_len = 0.
        count = 0
        while len_tip - prev_len >0.:
            #print count
            count+= 1
            prev_len += segment_length
            edget = '+' if cid == _root_id else '<'
            cid = g.add_child(cid, edge_type=edget, label='S', base_length=prev_len+parent_base,
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
        np.random.seed(seed)

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
                    print 'shift ', shift, i
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
        _axis = list(algo.descendants(g,vid,RestrictedTo='SameAxis'))
        n = len(_axis)
        print 'AXIS ', _axis
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

        for index, (ramif, position) in enumerate(shuffled_axis):
            if ramif :

                aid = _axis[index+1]
                nid = g.node(aid)
                anchors.append(nid)

    if SUPERIOR_ORDER:
        # Update the ramification of order 1
        for v in g:
            if g.edge_type(v) == '+':
                print 'ORDER: ', v, g.order(v)
                update_randomized_delayed_axis(v)

        print 'ANCHORS ', anchors

        while anchors:   # while they are branching point left
            nid = anchors.pop(0)  # take next branching point
            if nid.order == 2:
                print nid, nid.position_index
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

    print 'branching_delay ', branching_delay
    print 'max_order', max(g.property('order').values())
    return g


def mtg_from_aqua_data(df, segment_length=1e-4):
    """ Added F. Bauget 2019-12-19
        reconstruct MTG from file in format used by aquaporin team
        maximum order is 2

        Author: C. Pradal
        Modified: F. Bauget


        Parameters
        ==========
            - df: pandas dataframe, 3 columns ['db','lr','order'], db and lr are length in m
            - segment_length: length of the vertices, default 1.e-4 in m

        Returns
        =======
            - g: MTG with the following properties set: edge_type, label, base_length, length

         the format is: 3 columns separated by tab
         * 1st col: "db" distance in m from base on the parent root where starts the lateral root
         * 2nd col: "lr" length in m of the corresponding lateral root
         * 3d col: "order" = 1 if parent root is the primary root, = 1-n if the parent root is a lateral root that starts at the node n on the parent root

    """

    g = MTG()
    rid = vid = g.add_component(g.root)

    rnid = g.node(rid)
    rnid.base_length = 0.
    rnid.length = segment_length
    rnid.label = 'S'
    rnid.order = 0

    # Primary root
    prev_len = 0.
    ramifs = {}  # variable (nb racine (1 first), index entry) [(index vertice node, lenght lat )]

    df_order = df[df.order == '1']  # array with 1st root
    length_base = df_order.index

    #    path = [1]
    count = 0
    for i in length_base:
        len_base = df_order.iloc[i].db
        code = '1'
        while len_base - prev_len > 0:
            # we add segment of segment_length till the next vertice => no lateral root yet so the edge_type is '<'
            prev_len += segment_length
            vid = g.add_child(vid, edge_type = '<', label = 'S', base_length = prev_len, length = segment_length,
                              order = 0, code = code)
        # Modification Decamber 2019by Fabrice: problem was that two following vertice may have the same base_length causing wrong g.property('position') recalculation in the conductance calculation
        # did not pass test_reconstruct_from_aqua_data in test_archi_data.py
        # vid = g.add_child(vid, edge_type='<', label='S', base_length=len_base, length=segment_length, order=0)
        # prev_len = len_base # the last added vertice may be > than the first node so we change len_base which is to big
        len_lateral = df_order.iloc[i].lr
        if len_lateral > 0.:
            count += 1
            p = tuple([1, count])  # 1: PR, count: countieme RL
            ramifs.setdefault(p, []).append((vid, len_lateral))  # randomly added, to sort it sorted(ramifs)

    ramifs = add_branching(g, df, ramifs = ramifs, Order = 1, segment_length = segment_length)

    ramifs = add_branching(g, df, ramifs = ramifs, Order = 3, segment_length = segment_length)
    return g

def add_branching(g, df, ramifs = None, Order = 0, segment_length = 1e-4):
    """ F. Bauget 2019-12-19
    add branching of a given order on the previous order
    linked to mtg_from_aqua_data
    Parameters:
        - g: MTG
        - df: pandas dataframe, 3 columns ['db','lr','order'], db and lr are length in m
        - ramifs: dict with a list (Order, nth lateral root), and a dict [vid, lr] vid is the vertice index on the parent root from which the lateral of length lr starts
        - Order: int the order of the new branching
        - segment_length: float length in m of the vertices

    Return:
        - new_ramifs: dict to used as the ramifs parameter for a new call of add_branching
    """
    new_ramifs = {}
    len_base = 0
    count = 0  # useless ?
    for path in ramifs:
        vid, lr = ramifs[path][
            0]  # vid is the vertice index on the parent root from which the lateral of length lr starts
        order = '-'.join(map(str, path))
        df_order = df[df.order == order]
        length_base = df_order.index
        code = order

        len_base = lr

        prev_len = 0.
        _root_id = vid
        parent_base = g.node(vid).base_length

        if df_order.empty:
            while len_base - prev_len > 0:
                prev_len += segment_length
                edge_type = '+' if _root_id == vid else '<'
                vid = g.add_child(vid, edge_type = edge_type, label = 'S', base_length = parent_base + prev_len,
                                  length = segment_length, order = Order, code = code)
        else:
            for i in length_base:
                len_base = df_order.db[i]
                edge_type = '+' if _root_id == vid else '<'
                while len_base - prev_len > 0:
                    prev_len += segment_length
                    vid = g.add_child(vid, edge_type = edge_type, label = 'S', base_length = parent_base + prev_len,
                                      length = segment_length, order = Order, code = code)
                    edge_type = '<'
                #Modification December 2019 by Fabrice: problem was that two following vertice may have the same base_length causing wrong g.property('position') recalculation in the conductance calculation
                # did not pass test_reconstruct_from_aqua_data in test_archi_data.py
                # vid = g.add_child(vid, edge_type='<', label='S', base_length=parent_base+len_base+segment_length, length=segment_length, order=1, code=code)
                # prev_len = len_base
                len_lateral = df_order.lr[i]
                if len_lateral > 0.:
                    count += 1
                    p = tuple([Order + 1, count])
                    new_ramifs.setdefault(p, []).append((vid, len_lateral))
    return new_ramifs