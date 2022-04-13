"""Analysis procedures for Hydroroot.

Define a set of methods to ease the analysis and simulation.
"""
import pandas
# from hydroroot.main import hydroroot
from openalea.mtg.traversal import pre_order2
from openalea.mtg.algo import orders


def nb_roots(g, l, root=1, dl=1e-4, max_order=None):
    """Compute the number of segment intercepted at a given length.

    Parameters
    ==========
        - g: MTG
        - l: length
        - root: the root vertex from which the tree traversal start
        - dl: length of MTG segments (use for g.property('mylength') calculation
        - max_order: maximum order of considering roots

    Returns
    =======
        - number of segment
    """
    length = {}

    if 'mylength' in g.property_names():
        length = g.property('mylength')
    else:
        for v in pre_order2(g, root):
            pid = g.parent(v)
            length[v] = length[pid] + dl if pid else dl
        g.properties()['mylength'] = length

    order = None
    if max_order is not None:
        if 'order' in g.property_names():
            order = g.property('order')
        else:
            order = orders(g, scale=g.max_scale())

    count = 0
    for v in g:

        if order and order.get(v) >= max_order:
            continue

        pid = g.parent(v)
        if pid and (length[pid] <= l <= length[v]):
            count += 1
    return count


def intercept(g, dists, dl=1e-4, max_order=None):
    # F. Bauget 2021-07-09 : added dl to arguments
    """
    Compute intercepts at given lengths from collet.

    Parameters
    ----------
        - g: (MTG)
        - dists: (list of Float) list distances from the collet
        - dl: (float) length of MTG segments
        - max_order: (int) maximum order of considering roots

    Returns
    -------
        - intercepts: list of number of intercepts according to distances in dists
    """


    intercepts = [nb_roots(g, x, dl = dl, max_order=max_order) for x in dists]
    return intercepts


def read_data(data):
    """Merge data and return a Dataframe."""
    names = ('relative_position', 'internode_length', 'LR_length', 'distance_to_tip')

    df = None
    for d in data:
        pd = pandas.read_csv(d, sep=';', header=1,
                             names=names)
        pd.sort('distance_to_tip', inplace=True)
        pd.LR_length.cumsum()
        pd['cum'] = pd.LR_length.cumsum()
        pd['length'] = pd.internode_length.cumsum()

        if df is None:
            df = pd
        else:
            df = df.append(pd)

        pd = df

    pd['cumsum'] = pd.LR_length.cumsum()
    pd.sort('distance_to_tip', inplace=True)

    return pd
