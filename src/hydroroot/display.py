import numpy as np
import matplotlib as mpl
from matplotlib import pyplot
from pylab import cm
from matplotlib.colors import Normalize, LogNorm

from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color

import openalea.plantgl.all as pgl

from .radius import discont_radius

def get_root_visitor(prune=None, factor = 1.0e4):
    """
    Turtle going through the architecture
    used in :func:`plot` and :func:`get_root_visitor_with_point`
    Root angles are dummy values to get better display in 3D

    :param prune: (float) - distance from base after witch the MTG is no longer read (Default value = None)
    :param factor: (float) - a factor apply to length and radius properties (Default value = 1.0e4)

    :For example: if the MTG vertices length is 0.1 mm (1.0e-4 m), then if we want to set each vertex
                  as a dot (let say a pixel) we have to multiply length by 1/1e-4 is useful for plot
    """
    # F. Bauget 2021-06-10 : added parameters factor_length and factor_radius because they were hard codded to 1.0e4
    # without any reason unless historically the segment_length was chosen = 1.0e-4
    # but may be a problem if we export in rsml because the unit may be wrong

    def root_visitor(g, v, turtle, prune=prune):
        """
        A visitor with different root angles according to their order for display purpose

        :param g: (MTG)
        :param v: (int) - node id
        :param turtle: - mtg.turtle.PglTurtle
        :param prune: (float) - distance from base after witch the MTG is no longer read (Default value = None)

        """
        mylength = {}
        if prune and ('mylength' in g.properties()):
            mylength = g.property('mylength')
        angles = [90,45]+[30]*5
        n = g.node(v)
        # radius = n.radius*1.e4
        radius = n.radius * factor
        order = int(n.order)
        # length = n.length*1.e4
        length = n.length * factor

        if prune:
            if mylength.get(v,0.)>prune:
                return

        if g.edge_type(v) == '+':
            angle = angles[order]
            turtle.down(angle)

        turtle.setId(v)
        turtle.setWidth(radius)
        for c in n.children():
            if c.edge_type() == '+':
                turtle.rollL(130)
        #turtle.setColor(order+1)
        turtle.F(length)

        # define the color property
        #n.color = random.random()
    return root_visitor

def plot_old(g, has_radius=False, r_base=1.e-4, r_tip=5e-5,
         visitor=None, prop_cmap='radius', cmap='jet',lognorm=False,
         prune=None):
    """
    :Deprecated:

    Create a scene from g

    :param g: MTG
    :param has_radius: boolean (Default value = False)
    :param r_base: float (Default value = 1.e-4)
    :param r_tip: float (Default value = 5e-5)
    :param visitor: Turtle going through the architecture see get_root_visitor (Default value = None)
    :param prop_cmap: property used for the color map (Default value = 'radius')
    :param cmap: matplotlib color map (Default value = 'jet')
    :param lognorm: True log scale normalization (Default value = False)
    :param prune: float (Default value = None)
    :returns: scene
    
    
    Exemple:

    >>> from openalea.plantgl.all import *
        >>> s = plot()
        >>> shapes = dict( (x.getId(), x.geometry) for x in s)
        >>> Viewer.display(s)
    """
    if visitor is None:
        visitor = get_root_visitor(prune=prune)

    r_base, r_tip = float(r_base), float(r_tip)

    if not has_radius:
        discont_radius(g,r_base=r_base, r_tip=r_tip)

    turtle = turt.PglTurtle()
    turtle.down(180)
    scene = turt.TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False)

    # Compute color from radius
    color.colormap(g,prop_cmap, cmap=cmap, lognorm=lognorm)

    # F. Bauget 2022-03-15: WIP python 2 to 3, got "AttributeError: 'Shape' object has no attribute 'getId'"
    # shapes = dict( (sh.getId(),sh) for sh in scene)
    shapes = scene.todict()
    colors = g.property('color')
    for vid in colors:
        if vid in shapes:
            # shapes[vid].appearance = pgl.Material(colors[vid])
            for sh in shapes[vid]:
                sh.appearance = pgl.Material(colors[vid])
    # scene = pgl.Scene(list(shapes.values()))
    scene = pgl.Scene([sh for shid in shapes.values() for sh in shid ])
    # F. Bauget 2022-03-15
    return scene

def plot(g = None, min=None, max=None, name=None, cmap = 'jet', **kwds):
    """
    Display the MTG g in 3D see :func:`mtg_scene`

    :param g: (MTG) - the mtg to plot (Default value = None)
    :param min: (float) - the minimum used for normalization (Default value = None)
    :param max: (float) - the maximum used for normalization (Default value = None)
    :param name: (string) - if not None save the plot to name (Default value = None)
    :param cmap: (string) - the color map name (https://matplotlib.org/stable/tutorials/colors/colormaps.html) (Default value = 'jet')
    :param kwds: see :func:`mtg_scene`
    """

    pgl.Viewer.display(mtg_scene(g, min = min, max = max, cmap = cmap, **kwds))
    if name is not None:
        pgl.Viewer.frameGL.saveImage(name)

def mtg_scene(g, has_radius=False, r_base=1.e-4, r_tip=5e-5,
         visitor=None, prop_cmap='radius', cmap='jet',lognorm=False,
         min = None, max = None, prune=None, factor = 1.0e4):
    """
    Build scene of a MTG to be displayed in 3D plantgl viewer (https://github.com/openalea/plantgl)

    :param g: (MTG)
    :param has_radius: (boolean) - if True use the \'radius\' property otherwise compute the radius (:func:`hydroroot.radius.discont_radius`) (Default value = False)
    :param r_base: (float) - radius (m) at the base to use if the radius property is computed (Default value = 1.e-4)
    :param r_tip: (float) - radius (m) at the tip to use if the radius property is computed (Default value = 5e-5)
    :param visitor: Turtle going through the architecture see :func:`get_root_visitor` (Default value = None)
    :param prop_cmap: (string) property used for the color map (Default value = 'radius')
    :param cmap: matplotlib color map (Default value = 'jet')
    :param lognorm: (boolean) - if None no normalization, if False linear normalization, if True log normalization (Default value = False)
    :param min: (float) - if not None, the minimum used for the colormap normalization (Default value = None)
    :param max: (float) - if not None, the maximum used for the colormap normalization (Default value = None)
    :param prune: (float) - distance from base after witch the MTG is no longer read (Default value = None)

    :returns: a plantgl Scene

    """
    # hydroroot.display.plot modified to use my_colormap where min-max can be imposed

    if visitor is None:
        visitor = get_root_visitor(prune=prune, factor = factor)

    r_base, r_tip = float(r_base), float(r_tip)

    if not has_radius:
        discont_radius(g,r_base=r_base, r_tip=r_tip)

    # Compute color from radius
    my_colormap(g,prop_cmap, cmap=cmap, lognorm=lognorm, min = min, max = max)

    turtle = turt.PglTurtle()
    turtle.down(180)
    scene = turt.TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False)

    # F. Bauget 2022-08-24 : lines below are duplicates with lines in turt.TurtleFrame
    # shapes = scene.todict()
    # colors = g.property('color')
    # for vid in colors:
    #     if vid in shapes:
    #         for sh in shapes[vid]:
    #             sh.appearance = pgl.Material(colors[vid])
    #
    # scene = pgl.Scene([sh for shid in shapes.values() for sh in shid ])

    return scene

def my_colormap_old(g, property_name, cmap='jet',lognorm=True):
    """

    :param g: 
    :param property_name: 
    :param cmap:  (Default value = 'jet')
    :param lognorm:  (Default value = True)

    """
    prop = g.property(property_name)
    keys = list(prop.keys())
    values = np.array(list(prop.values()))
    #m, M = int(values.min()), int(values.max())
    _cmap = cm.get_cmap(cmap)
    norm = Normalize() if not lognorm else LogNorm()
    values = norm(values)
    #my_colorbar(values, _cmap, norm)

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(list(zip(keys,colors)))

def my_colormap(g, property_name, cmap='jet',lognorm = False, min = None, max = None):
    """
    Compute the property \`color\` based on a given property and a colormap.

    :param g: (MTG)
    :param property_name: (string) - the property to display with  the heatmap
    :param cmap: (string) - the color map name (https://matplotlib.org/stable/tutorials/colors/colormaps.html) (Default value = 'jet')
    :param lognorm: if True Normalize to the 0-1 range on a log scale, if False on a linear scale (Default value = False)
    :param min: (float) - if not None, the minimum used for normalization (Default value = None)
    :param max: (float) - if not None, the maximum used for normalization (Default value = None)

    :returns: - g (MTG) with the property 'color' set


    .. note::
        This is the function openalea.mtg.plantframe.color.colormap modified to add the possibility to fix the min-max
        for normalization and so the possibility to compare two MTG with different min-max but the same colormap
    """
    # F. Bauget 2022-07-26: the previous function my_colormap_old modified with the addition of min and max argument
    prop = g.property(property_name)
    keys = list(prop.keys())
    values = np.array(list(prop.values()))

    _cmap = color.get_cmap(cmap)
    if lognorm == False: # explecitly == False
        norm = Normalize(vmin = min, vmax = max)
        values = norm(values)
    elif lognorm == True: # explecitly == True
        norm = LogNorm(vmin = min, vmax = max)
        values = norm(values)
    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(list(zip(keys,colors)))
    return g

def my_colorbar_old(values, cmap, norm):
    """
    display a color scale

    :param values: 
    :param cmap: 
    :param norm: 

    """
    fig = pyplot.figure(figsize=(8,3))
    ax = fig.add_axes([0.05, 0.65, 0.9, 0.15])
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap, norm=norm, values=values)

def property_scale_bar(g, property_name, cmap='jet',lognorm = False, vmin = None, vmax = None, format=None):
    """
    display a color scale bar based on the property_name values

    :param g: (MTG)
    :param property_name: (string) - the property name
    :param cmap: (string) - the color map name (Default value = 'jet')
    :param lognorm: (boolean) - if True Normalize to the 0-1 range on a log scale (Default value = False)
    :param min: (float) - if not None the scale minimum (Default value = None)
    :param max: (float) - if not None the scale maximum (Default value = None)
    :param format: (string) - format for tick labes, e.i. '%.2f' (Default value = None)
    """
    prop = g.property(property_name)
    values = np.array(list(prop.values()))

    if vmin is None:
        vmin = values.min()
    if vmax is None:
        vmax = values.max()

    if lognorm == False:
        norm = Normalize(vmin = vmin, vmax = vmax)
        # values = norm(values)
    else:
        norm = LogNorm(vmin = vmin, vmax = vmax)
        # values = norm(values)

    fig, ax = pyplot.subplots(figsize = (6, 1))
    fig.subplots_adjust(bottom = 0.5)
    cb1 = mpl.colorbar.Colorbar(ax, cmap = cmap, norm = norm, orientation = 'horizontal', format = format)
    cb1.set_label(property_name)
    fig.show()

def get_root_visitor_with_point(prune=None, factor = 1.0e4):
    """Get 3D position of root segment by using turtle in :func:`get_root_visitor`
    Create a new property position3d = [x,y,z] coordinate needed for example to export MTG in RSML format see
    :func:`hydroroot.hydro_io.export_mtg_to_rsml`

    :param prune: (float)  - distance from base after witch the MTG is no longer read (Default value = None)
    :param factor: (float) a factor apply to property length (Default value = 1.0e4)

    .. note::
        the 3D coordinate are calculated from a virtual 3D representation with imaginary angles. Therefore, the resulting
        3D coordinate are dummy values.

    """
    visitor = get_root_visitor(prune=prune, factor = factor)

    def root_visitor3D(g, v, turtle, prune=prune):
        """
        Add position3d to the visitor from :func:`get_root_visitor`

        :param g: (MTG)
        :param v: (int) - node id
        :param turtle: - mtg.turtle.PglTurtle
        :param prune: (float) - distance from base after witch the MTG is no longer read (Default value = None)

        """
   
        visitor(g, v, turtle, prune=prune)
        n = g.node(v)
        n.position3d = tuple(turtle.getPosition())

    return root_visitor3D
