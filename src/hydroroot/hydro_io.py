"""
author: F. Bauget
date: 2019-12-20

Functions related to io process in hydroroot
"""

from openalea.mtg.algo import *
from openalea.mtg.traversal import *
import pandas as pd


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

    v_base = g.component_roots_at_scale_iter(g.root, scale = g.max_scale()).next()

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
                if g.edge_type(v2) == '+' and parent <> parent2:
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
