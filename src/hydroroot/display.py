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

def get_root_visitor(prune=None):
    def root_visitor(g, v, turtle, prune=prune):
        mylength = {}
        if prune and ('mylength' in g.properties()):
            mylength = g.property('mylength')
        angles = [90,45]+[30]*5
        n = g.node(v)
        radius = n.radius*1.e4
        order = int(n.order)
        length = n.length*1.e4

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

    shapes = dict( (sh.getId(),sh) for sh in scene)

    colors = g.property('color')
    for vid in colors:
        if vid in shapes:
            shapes[vid].appearance = pgl.Material(colors[vid])
    scene = pgl.Scene(shapes.values())
    return scene


def my_colormap(g, property_name, cmap='jet',lognorm=True):
    prop = g.property(property_name)
    keys = prop.keys()
    values = np.array(prop.values())
    #m, M = int(values.min()), int(values.max())
    _cmap = cm.get_cmap(cmap)
    norm = Normalize() if not lognorm else LogNorm()
    values = norm(values)
    #my_colorbar(values, _cmap, norm)

    colors = (_cmap(values)[:,0:3])*255
    colors = np.array(colors,dtype=np.int).tolist()

    g.properties()['color'] = dict(zip(keys,colors))


def my_colorbar(values, cmap, norm):
    fig = pyplot.figure(figsize=(8,3))
    ax = fig.add_axes([0.05, 0.65, 0.9, 0.15])
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap, norm=norm, values=values)


def get_root_visitor_with_point(prune=None):
    visitor = get_root_visitor(prune=prune)

    def root_visitor3D(g, v, turtle, prune=prune):
   
        visitor(g, v, turtle, prune=prune)
        n = g.node(v)
        n.position3d = tuple(turtle.getPosition())

    return root_visitor3D
