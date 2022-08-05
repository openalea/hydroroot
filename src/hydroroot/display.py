from math import pi

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot
from pylab import cm, colorbar
from pylab import plot as pylab_plot
from matplotlib.colors import Normalize, LogNorm

from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color

import openalea.plantgl.all as pgl

from .radius import discont_radius

def get_root_visitor(prune=None, factor = 1.0e4):
    # F. Bauget 2021-06-10 : added parameters factor_length and factor_radius because they were hard codded to 1.0e4
    # without any reason unless historically the segment_length was chosen = 1.0e-4
    # but may be a problem if we export in rsml because the unit may be wrong
    """
    Turtle going through the architecture
    used in plot and get_root_visitor_with_point

    :parameters:
        - `prune` (float) - distance from base after witch the MTG is no longer read
        - `factor` - a factor apply to length and radius properties

    :For example: if the MTG vertices length is 0.1 mm (1.0e-4 m), then if we want to set each vertex
                  as a dot (let say a pixel) we have to multiply length by 1/1e-4 is useful for plot
    """
    def root_visitor(g, v, turtle, prune=prune):
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

def plot(g, has_radius=False, r_base=1.e-4, r_tip=5e-5,
         visitor=None, prop_cmap='radius', cmap='jet',lognorm=False,
         prune=None):
    """
    Deprecated look at
    Exemple: mtg_scene

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

def mtg_scene(g, has_radius=False, r_base=1.e-4, r_tip=5e-5,
         visitor=None, prop_cmap='radius', cmap='jet',lognorm=False,
         min = None, max = None, prune=None):
    """
    hydroroot.display.plot modified to use my_colormap where min-max can be imposed

    :Parameters:
		- g: (MTG)
		- has_radius: (boolean) True use of g.property('radius'), False radii are calculated
		- r_base: (float) radius of the base (m) used for radius calculation if has_radius = True
		- r_tip: (float) radius of the tip (m) used for radius calculation if has_radius = True
		- visitor: Turtle going through the architecture see get_root_visitor
		- prop_cmap: property used for the color map
		- cmap: matplotlib color map
		- lognorm: True log scale normalization, False normal normalization used for color map
		- prune: (float) - distance from base after witch the MTG is no longer read

    :return: 
		- scene a plantgl Scene
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
    my_colormap(g,prop_cmap, cmap=cmap, lognorm=lognorm, min = min, max = max)

    shapes = scene.todict()
    colors = g.property('color')
    for vid in colors:
        if vid in shapes:
            for sh in shapes[vid]:
                sh.appearance = pgl.Material(colors[vid])

    scene = pgl.Scene([sh for shid in shapes.values() for sh in shid ])

    return scene

def my_colormap_old(g, property_name, cmap='jet',lognorm=True):
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

def my_colormap(g, property_name, cmap='jet',lognorm=True, min = None, max = None):
    # F. Bauget 2022-07-26: the previous function my_colormap_old modified with the addition of min and max argument
    """
    Compute a color property based on a given property and a colormap.
    openalea.mtg.plantframe.color.colormap modified to add the possibility to fix the min-max for normalization
    and so the possibility to compare two MTG with different min-max but the same colormap
    :Parameters:
    	- g: (MTG)
    	- property_name: (string) - the property to display with  the heatmap
    	- cmap: (string) - the heatmap name (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
    	- lognorm: (boolean) - False linear normalization, log normalization otherwise
    	- min: (float) - the minimum used for normalization
    	- max: (float) - the maximum used for normalization
    :Returns:
        - g (MTG) with the property 'color' set
    """

    prop = g.property(property_name)
    keys = list(prop.keys())
    values = np.array(list(prop.values()))

    _cmap = color.get_cmap(cmap)
    norm = Normalize(vmin = min, vmax = max) if not lognorm else LogNorm(vmin = min, vmax = max)
    values = norm(values)

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(list(zip(keys,colors)))
    return g

def my_colorbar(values, cmap, norm):
    fig = pyplot.figure(figsize=(8,3))
    ax = fig.add_axes([0.05, 0.65, 0.9, 0.15])
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap, norm=norm, values=values)


def get_root_visitor_with_point(prune=None, factor = 1.0e4):
    """
    Get 3D position of root segment by using turtle in get_root_visitor
    Create a new property position3d = [x,y,z] coordinate

    :parameters:
        - `prune` (float) - distance from base after witch the MTG is no longer read
        - `factor` - a factor apply to property length (see get_root_visitor)

    Remark: the 3D coordinate are calculated from an virtual 3D representation with imaginary angles
    """
    visitor = get_root_visitor(prune=prune, factor = factor)

    def root_visitor3D(g, v, turtle, prune=prune):
   
        visitor(g, v, turtle, prune=prune)
        n = g.node(v)
        n.position3d = tuple(turtle.getPosition())

    return root_visitor3D
