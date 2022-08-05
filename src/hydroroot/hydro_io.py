"""
author: F. Bauget
date: 2019-12-20

Functions related to io process in hydroroot
"""
import pandas as pd
import numpy as np

from openalea.mtg import MTG
from openalea.mtg.algo import axis
from openalea.mtg.traversal import *

from rsml import continuous, io

from hydroroot import display



def export_mtg_to_aqua_file(g, filename = "out.csv"):
    """
    Export a MTG architecture in a csv file into format used by aquaporin team
    g: MTG
    filename: the name of the output file

          the format is: 3 columns separated by tab
         * 1st col: "distance_from_base_(mm)" distance in mm from base on the parent root where starts the lateral root
         * 2nd col: "lateral_root_length_(mm)" length in mm of the corresponding lateral root
         * 3d col: "order" = 1 if parent root is the primary root, = 1-n if the parent root is a lateral root that
                            starts at the node n on the parent root

    It uses only the mtg properties 'position', 'order' and 'edge_type' because they are the only ones saved in
    simulated architecture
    At this stage (2019-12-20) only up to the 2d order
    """

    results = {'distance_from_base_(mm)': [], 'lateral_root_length_(mm)': [], 'order': []}

    v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))

    count = 0
    for v in pre_order2(g, 1):
        parent = g.parent(v)
        if g.edge_type(v) == '+' and g.order(v) == 1:
            # position is the length from tip so the total racine length is the 1st position on the racine
            racine_length = g.property('position')[min(axis(g, parent))] * 1e3  # unit change: m to mm
            results['distance_from_base_(mm)'].append(racine_length - g.property('position')[parent] * 1e3)

            LR_length = g.property('position')[v] * 1e3
            results['lateral_root_length_(mm)'].append(LR_length)
            results['order'].append(1)

    results['distance_from_base_(mm)'].append(racine_length)
    results['lateral_root_length_(mm)'].append(0)
    results['order'].append(1)

    count = 0
    count2 = 0
    for v in pre_order2(g, 1):
        parent = g.parent(v)
        if g.edge_type(v) == '+' and g.order(v) == 1:
            count2 = 0
            count += 1
            for v2 in pre_order2(g, v):
                parent2 = g.parent(v2)
                if g.edge_type(v2) == '+' and parent != parent2:
                    racine_length = g.property('position')[min(axis(g, parent2))] * 1e3
                    results['distance_from_base_(mm)'].append(racine_length - g.property('position')[parent2] * 1e3)

                    LR_length = g.property('position')[v2] * 1e3
                    results['lateral_root_length_(mm)'].append(LR_length)
                    results['order'].append('-'.join((str(1), str(count))))

                    count2 = count

            if count == count2:
                results['distance_from_base_(mm)'].append(racine_length)
                results['lateral_root_length_(mm)'].append(0)
                results['order'].append('-'.join((str(1), str(count))))

    df = pd.DataFrame(results, columns = ['distance_from_base_(mm)', 'lateral_root_length_(mm)', 'order'])
    df.to_csv(filename, sep = '\t', index = False)

def export_mtg_to_rsml(g_discrete, filename = None, segment_length = 1.0e-4):
    """
    Export a discrete MTG architecture into a rsml file or return a continuous MTG according to the parameters (see below)
    only the geometry is exported no other properties

    :Parameters:
        - `g_c` (MTG) - the discrete MTG to export
        - `filename` (string) - the name of the output file, if None return a conitnuous MTG
        - `segment_length` (float) - length of the MTG segments in meter

    Remark: g_discrete from hydroroot has length and radius in (m)

    - Does not overwrite the MTG in input
    - 1st: use a turtle to get position in 3D
            - use functions from hydroroot.display
            - in display.get_root_visitor_with_point() g property "position3d" is created
    - 2d: insert scales
            - MTG from hydroroot has only segment scale, the finest
            - add the axes scale
            - add the plant scale
            - end up with: plant scale = 1, axes scale = 2, segments scale = 3
    - 3d: convert the discrete MTG to continuous
            - MTG without the finest scale, and with polylines to discribe axes
    - 4th: write to the file
    """

    g = g_discrete.copy()

    # HydroRoot MTG lengths are in meter and RSML are in pixel => 1 segment is a pixel
    # the resolution property in RSML metadata is the conversion factor from pixel to unit then here
    # resolution * pixel = segment_length => resolution = segment_length (m)
    factor = 1.0 / segment_length

    # Compute 3D polylines: new property 'position3d'
    visitor = display.get_root_visitor_with_point(factor = factor)
    scene = display.mtg_scene(g, visitor = visitor)  # improved display.plot and changed its name because it was not a plot but a scene

    # Scale insertion: axes
    def quotient_axis(v):
        if g.edge_type(v) == '+':
            return True
        elif g.parent(v) is None:
            return True
        else:
            return False
    g = g.insert_scale(inf_scale = 1, partition = quotient_axis, default_label = 'A')

    # scale insertion: plant
    root_id = next(g.component_roots_iter(g.root))
    g = g.insert_scale(inf_scale = 1, partition = lambda v: v == root_id, default_label = 'P')

    # Transform and write to rsml file
    g = continuous.discrete_to_continuous(g, position = 'position3d')

    # set some metadata
    g.graph_properties()['metadata']['unit'] = 'm'
    g.graph_properties()['metadata']['resolution'] = segment_length # pixel -> m
    g.graph_properties()['metadata']['software'] = 'HydroRoot'

    if filename is not None:
        io.mtg2rsml(g, filename)
        return
    else:
        return g

def import_rsml_to_discrete_mtg(g_c, segment_length = 1.0e-4, resolution = 1.0e-4):
    # F. Bauget 2020-03-18 : RSML continuous from rsml2mtg()  to hydroroot disctrete copied from rsml
    # don't use parent node because rsml from other places don't have them but only coordinates of polylines
    """
    Convert MTG from continuous to discrete usable in HydroRoot

    Does the reverse of `discrete_to_continuous`:
      - Add a sequence of segments to all axes from their `geometry` attribute


    Based only on the coordinates of the polylines:
        - calculation of the euclidian distance between points of polylines then add the segment needed
        - to place a lateral on its parent axe, we calcul the distance between the 1st point of the lateral
          and the vertices of the parent axe and take the shortest distance to select the branshing vertex on the axe


    :Parameters:
        - `g_c` (MTG) - the continuous MTG to convert
        - `segment_length` (Float) - the segment length in meter (m)
        - `resolution` (float) - the resolution of the polylines coordinates, in polylines unit per meter (unit/m)
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
    rnid.edge_type = '<'


    axe_segments = {}

    for axe in continuous.toporder(g_c, g_c.max_scale()):
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
        _length = np.sqrt(np.sum(vec, axis = 1)) * resolution


        if _order == 0:
            # primary root
            axe_segments[(axe, 0)] = seg # the first
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

        axe_segments[(axe, 1 )] = seg

        # create the other segments
        for i, l in enumerate(_length[1:]):
            if l > 0.0:
                n = int(l/segment_length)
                if n == 0: n = 1
            for j in range(int(n)):
                seg = g.add_child(seg, edge_type = '<', label = 'S', length = segment_length, order = _order)
            axe_segments[(axe, i + 2)] = seg


    return g
