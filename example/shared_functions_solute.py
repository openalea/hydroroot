import numpy as np
import pandas as pd
import math


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize, LogNorm

from openalea.mtg import traversal, algo
from openalea.mtg.algo import axis
from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color
from openalea.plantgl.all import Viewer, Scene, Material
from openalea.mtg.io import write_mtg

from hydroroot.display import get_root_visitor
from hydroroot.radius import discont_radius



def to_mtg_text_file(g, filename = 'example.mtg'):
    # Export all the properties defined in `g`.
    # We consider that all the properties are real numbers.

    properties = []
    for p in g.property_names():
        if p in ['edge_type', 'code', 'label']:
            properties.append((p, 'TXT'))
        else:
            properties.append((p, 'REAL'))
    mtg_lines = write_mtg(g, properties)

    # Write the result into a file example.mtg
    f = open(filename, 'w')
    f.write(mtg_lines)
    f.close()

def export_mtg_to_aqua_data(g):
    """
    Export a MTG architecture in a csv file into format used by aquaporin team

    :param g: (MTG) - the root architecture

    :return: df: (pandas DataFrame) - 3 columns

          the format is: 3 columns separated by tab
         * 1st col: "distance_from_base_(mm)" distance in mm from base on the parent root where starts the lateral root
         * 2nd col: "lateral_root_length_(mm)" length in mm of the corresponding lateral root
         * 3d col: "order" = 1 if parent root is the primary root, = 1-n if the parent root is a lateral root that
                            starts at the node n on the parent root

    It uses only the mtg properties 'position', 'order' and 'edge_type' because they are the only ones saved in
    simulated architecture
    At this stage (2019-12-20) only up to the 2d order
    """
    # Activate(g)

    results = {'distance_from_base_(mm)': [], 'lateral_root_length_(mm)': [], 'order': []}

    v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))

    count = 0
    for v in traversal.pre_order2(g, 1):
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
    for v in traversal.pre_order2(g, 1):
        parent = g.parent(v)
        if g.edge_type(v) == '+' and g.order(v) == 1:
            count2 = 0
            count += 1
            for v2 in traversal.pre_order2(g, v):
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
    return df

def plot(g = None, filename=None, min=None, max=None, name=None, cmap = 'jet', **kwds):
    """
    plot a MTG in 3D given in argument or from a file

    :param g: (MTG) - the mtg to plot
    :param filename: (string) - if not None plot the architecture given in filename in the Aqua team format
    :param min: (float) - the minimum used for normalization
    :param max: (float) - the maximum used for normalization
    :param name: (string) - if not None save the plot to name
    :param cmap:  (string) - the heatmap name (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
    :param kwds: has_radius=False, r_base=1.e-4, r_tip=5e-5,
             visitor=None, prop_cmap='radius', cmap='jet',lognorm=False,
             prune=None
    :return:
    """

    if filename is not None:
        df = read_archi_data(filename)
        g, primary_length, _length, surface, _seed = root_builder(df = df)

    Viewer.display(my_mtg_scene(g, min = min, max = max, cmap = cmap, **kwds))
    # Viewer.display(mtg_scene(g, cmap = cmap, **kwds))
    if name is not None:
            Viewer.frameGL.saveImage(name)

def my_plot_with_bar(g, prop, min=None, max=None, lognorm = False, cmap = 'jet', format=None, name=None, **kwds):

    """
    plot a MTG in 3D given in argument or from a file and the heatmap scale bar

    :param g: (MTG) - the mtg to plot
    :param prop: (string) - the property to plot
    :param min: (float) - the minimum used for normalization
    :param max: (float) - the maximum used for normalization
    :param lognorm: (boolean) -  True log scale normalization, False normal normalization used for color map
    :param cmap: (string) - the heatmap name (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
    :param format: (str or Formatter) - the labels format
    :param name: (string) - if not None save the plot to name
    :param kwds: has_radius=False, r_base=1.e-4, r_tip=5e-5,
             visitor=None, prop_cmap='radius', cmap='jet',lognorm=False,
             prune=None
    :return:
    """

    values = np.array(list(g.property(prop).values()))
    plot(g, min = min, max = max, prop_cmap = prop, lognorm = lognorm, cmap = cmap, name = name, **kwds)
    fig, ax = plt.subplots(figsize = (6, 1))
    fig.subplots_adjust(bottom = 0.5)

    if min is None: min = values.min()
    if max is None: max = values.max()
    if lognorm:
        norm = mpl.colors.LogNorm(vmin = min, vmax = max)
    else:
        norm = mpl.colors.Normalize(vmin = min, vmax = max)
    _cmap = color.get_cmap(cmap)
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap = _cmap, norm = norm, orientation = 'horizontal', format=format)
                                    # ,ticks=np.arange(0,0.13,0.01))
    cb1.set_ticks([min, 0.5*(min+max), max])
    cb1.set_label(prop)
    fig.show()
    if name is not None:
        fig.savefig('bar-'+name)

def my_plot_property(g, prop, max_order = 10):
    """

    :param g: g (MTG)
    :param prop: (string) - property to plot
    :param max_order: (int) - the maximum order to plot
    :return:
    :Example:
        >>> g.plot_property('length')
    """

    import matplotlib
    import matplotlib.pyplot
    import numpy as np
    props = g.property(prop)
    if max_order <= 10:
        pylab_colors = list(matplotlib.colors.TABLEAU_COLORS.keys())
    else:
        pylab_colors = list(matplotlib.colors.cnames.keys())
    color = {}
    orders = algo.orders(g)
    for k in props:
        color[k]=pylab_colors[orders[k]]
    #color = {k:pylab_colors[orders[k]] for k in props} # not in Py2.6

    heights = algo.heights(g)
    h = np.array([heights[v] for v in props])
    _prop = np.array(list(props.values()))
    for v in props:
        matplotlib.pyplot.plot(heights[v], props[v], 'o', color = color[v])

def my_toporder(g, scale):
    """ Return the list of `g` vertices at scale `scale` in topological order """
    axes = []
    list(map(axes.extend,(traversal.pre_order2(g,vid)
                          for vid in g.vertices(scale=scale)
                          if not g.parent(vid))))
    return axes

def viscosity_peg_old(Cpeg = 0.0, unit_factor = 1.0, T = 298):
    """
    Dynamic viscosity, mu, calculation in mPa.s according to the PEG-8000 concentration for a temperature T
    equation 5 from Gonzllez-Tello J. Chem. Eng. Data 1994,39, 611-614

    mu = 1.3311*exp((120.2+1143*w)/(T-172))
    mu in mPa.s
    T in Kelvin
    w in g/g water

    for w e [0 ; 0.1] : linear relation between points 1 and 1.3311*exp((120.2)/(T-172))
    This to get 1 mPa.s when w=0 (~ water viscosity)

    parameters:
    - Cpeg: Float (0), PEG concentration
    - unit_factor: Float (1), factor used to pass the concentration in g/g water unit
    - T: Float (293), temperature in Kelvin

    return:
    - mu: Float, the dynamic viscosity in mPa.s
    """

    if Cpeg <= 1.0e-20: Cpeg = 0.0
    # if Cpeg > 150e-9: Cpeg=150e-9
    w = Cpeg * unit_factor
    if w > 1.: w=1.0
    if w < 0.1:
        a=1.3311*math.exp((120.2+1143*0.1)/(T-172))
        mu = 1.0 + (a - 1.0) / 0.1 * w
    else:
        mu = 1.3311*math.exp((120.2+1143*w)/(T-172))

    return mu

def viscosity_peg6000(Cpeg = 0.0, unit_factor = 1.0):
    """
    Dynamic viscosity, mu, calculation in mPa.s according to the PEG concentration for a temperature of 293 K
    based on data from Mei et al., J. Chem. Eng. Data 1995, 40, 1168-1171
    eq. 1: mu = mu0 * exp(17.284*x-20.887*x^2) here the water viscosity mu0 ~ 1 mPa.s

    parameters:
    - Cpeg: Float (0), PEG concentration
    - unit_factor: Float (1), factor used to pass the concentration in %wt of water

    return:
    - mu: Float, the dynamic viscosity in mPa.s
    """

    if Cpeg < 1.0e-20: Cpeg = 0.0
    if Cpeg > 150e-9: Cpeg=150e-9
    x = Cpeg * unit_factor
    mu = math.exp(17.284*x-20.887*x**2)

    return mu

def viscosity_peg1000(Cpeg = 0.0, unit_factor = 1.0):
    """yyy
    Dynamic viscosity, mu, calculation in mPa.s according to the PEG concentration for a temperature of 293 K
    based on data from Mei et al., J. Chem. Eng. Data 1995, 40, 1168-1171
    A data point mu=1 mPa.s at Cpeg = 0 has been added to get the following law

    parameters:
    - Cpeg: Float (0), PEG concentration
    - scale_factor: Float (1), factor used to pass the concentration mol/m3

    return:
    - mu: Float, the dynamic viscosity in mPa.s
    """

    if Cpeg < 1.0e-20: Cpeg = 0.0
    if Cpeg > 150e-9: Cpeg=150e-9
    mu = 0.13 + 0.87 * math.exp(Cpeg * unit_factor / 153.0)

    return mu

def osmotic_p_peg6000(Cpeg = 0.0, unit_factor = 1.0):
    """
    osmotic pressure calculation of the PEG 6000 according to Steuder1981
    Y = 2.0e-5 X**1.94
    [Y] = MPa
    [X] = g/Kg of water or g/L

    parameters:
    - Cpeg: Float (0), PEG concentration
    - unit_factor: Float (1), factor used to pass the concentration in g/Kg of water

    return:
    - osmotic_p: Float, the osmotic pressure (MPa)
    """
    if Cpeg < 1.0e-20: Cpeg = 1.0e-20
    osmotic_p = - 1.0e-5 * (Cpeg * unit_factor)**2.0

    # osmotic_p = - 8.31 * 293 * Cpeg * 1e3

    return osmotic_p

def derivative_osmotic_p_peg6000(Cpeg = 0.0, unit_factor = 1.0):
    """
    the 1st derivative according to X of the osmotic pressure calculation of the PEG 6000 according to Steuder1981
    dY/dX = 2.0e-5 * 1.94 X**0.94
    [X] = g/Kg of water or g/L

    parameters:
    - Cpeg: Float (0), PEG concentration
    - unit_factor: Float (1), factor used to pass the concentration in g/Kg of water

    return:
    - derivative: Float, the osmotic pressure derivative according to Cpeg
    """
    if Cpeg < 1.0e-20: Cpeg = 1.0e-20
    derivative = - 2.0e-5 * (Cpeg * unit_factor) * unit_factor # to get MPa/(unit Cpeg)

    # derivative = - 8.31 * 293 * 1e3

    return derivative

def cbase_vs_jv(jv):
    """
    calculation of the exudate solute concentration according to the water flux
    the equation used is an exponential fit done on data of experiment 2021-08-24-FB-Exp06

    parameters:
    - jv: Float, the water flow in microL/s

    return:
    - c: Float, the concentration in mol/m3
    """
    c = 17+277*np.exp(-jv/0.00807)

    return c


