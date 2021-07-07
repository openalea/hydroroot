import numpy as np

from openalea.mtg import *
from hydroroot import display
from path import Path
from rsml.continuous import toporder

data = Path('./data/200703-YBFB-Col-2.mtg')
# gini = MTG(data)
g = MTG(data)

# 2 Compute 3D polylines
visitor = display.get_root_visitor_with_point()
scene = display.plot(g, visitor=visitor)

# Position is position3D

# 1 Scale insertion
def quotient_axis(v):
    if g.edge_type(v) == '+':
        return True
    elif g.parent(v) is None:
        return True
    else:
        return False

g = g.insert_scale(inf_scale=1, partition=quotient_axis, default_label='A')

root_id = next(g.component_roots_iter(g.root))
g = g.insert_scale(inf_scale=1, partition=lambda v: v==root_id, default_label='P')

#3 Transform to rsml 
from rsml import continuous, io, plot

gcopy = g.copy()

g = continuous.discrete_to_continuous(g, position='position3d')
io.mtg2rsml(g, 'myfile.rsml')

# #plot mtg from rsml file
# g2 = io.rsml2mtg('myfile.rsml')

def my_continuous_to_discrete(g_c, segment_length = 1.0e-4, resolution = 1.0e-4):
    # F. Bauget 2020-03-18 : testing RSML continuous from rsml2mtg()  to hydroroot disctrete copied from rsml
    # don't use parent node because rsml from other places don't have them but only coordinates of polylines
    """ Convert mtg `g` from continuous to discrete

    Does the reverse of `discrete_to_continuous`:
      - Add a sequence of segments to all axes from their `geometry` attribute


    Based only on the coordinate of the polylines:
        - calculation of the euclidian distance between points of polylines then add the segment needed
        - to place a lateral on its parent axe, we calcul the distance between the 1st point of the lateral
          and the vertices of the parent axe and take the shortest distance to select the branshing vertex on the axe


    :Parameters:
        - `g_c` (MTG) - the continuous MTG to convert
        - `segment_length` (Float) - the segment length in meter
        - `resolution` (float) - the resolution of the polylines coordinates, usually 1 pixels
    """


    geometry = g_c.property('geometry')

    _order = 0 # _order = 0 => 1st axe <=> primary root

    g = MTG()
    rid = seg = g.add_component(g.root)

    rnid = g.node(rid)
    rnid.base_length = 0.
    rnid.length = segment_length
    rnid.label = 'S'
    rnid.order = _order


    axe_segments = {}

    for axe in toporder(g_c, g_c.max_scale()):
        # get segment on parent branch
        p_axe = g_c.parent(axe)
        v = axe
        _order = 0
        while g_c.parent(v) is not None:
            v = g_c.parent(v)
            _order += 1


        if p_axe is not None:
            # parent axe
            pos_axe = np.array(geometry[p_axe])
            # vec_axe = np.diff(pos_axe,axis=0)**2
            # _length_axe = np.sqrt(np.sum(vec_axe, axis = 1))*resolution

        pos = np.array(geometry[axe])
        vec = np.diff(pos,axis=0)**2
        _length = np.sqrt(np.sum(vec, axis = 1))*resolution


        if _order == 0:
            # primary root
            if _length[0] > 0.0:
                n = int(_length[0]/segment_length)
                if n == 0: n = 1
            for j in range(n):
                seg = g.add_child(seg, edge_type='<', label = 'S', length = segment_length, order = _order)

        else:
            min = 1.0e10
            for i, p in enumerate(pos_axe):
                # search for the smallest distance between each nodes of p_axe and the 1st node of the child
                distance = np.sqrt(np.sum((p-pos[0])**2.0))
                if distance < min:
                    min = distance
                elif i > 1:
                    break

            seg = axe_segments[(p_axe, i - 1)] # branching vertex on parent axe

            seg = g.add_child(seg, edge_type='+', label = 'S', length = segment_length, order = _order)
            if _length[0] > 0.0:
                n = int(_length[0]/segment_length)
                if n == 0: n = 1
            for j in range(n-1):
                seg = g.add_child(seg, edge_type='<', label = 'S', length = segment_length, order = _order)

        shift = 0
        axe_segments[(axe, shift)] = seg

        # create the other segments
        for i, l in enumerate(_length[(1 + shift):]):
            if l > 0.0:
                n = int(l/segment_length)
                if n == 0: n = 1
            for j in range(int(n)):
                seg = g.add_child(seg, edge_type = '<', label = 'S', length = segment_length, order = _order)
            axe_segments[(axe, i + 1)] = seg


    return g
