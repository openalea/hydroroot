import random
import numpy as np

from math import pi
from matplotlib import pyplot, mpl
from pylab import cm, colorbar
from pylab import plot as pylab_plot
from matplotlib.colors import Normalize, LogNorm

from openalea.mtg import turtle as turt
import openalea.plantgl.all as pgl

from hydroroot.flux import *
from hydroroot import radius, conductance, markov



def surface(g, length=1e-4):
    surf = 0
    radius = g.property('radius')
    for vid in g.vertices():
        if radius.has_key(vid):
            surf += 2 * pi * radius[vid] * length
    return surf

def compute_flux(g, n=300, psi_e=300000., psi_base=101325., Jv=1e-10, k0=0.5, length=1e-4):
    k0 = float(k0)
    radius.discont_radius(g, r_base=1.e-4, r_tip=5.e-5)
    surf = surface(g)
    print 'surf',surf
    k = conductance.compute_k(g, k0,length)
    K = conductance.compute_K(g,length=length)

    assert all(v>0 for v in K.values()),K

    g = flux(g, k, K, Jv, psi_e, psi_base)

    J_out = g.property('J_out')
    assert all(v>0 for v in J_out.values()), J_out.values()
    return g

def test_linear(n=300, psi_e=300000., psi_base=101325., Jv=1e-10, k0=0.5, length=1e-4):
    """ Test flux and water potential computation on a linear root.

    Units :
    psi_e & psi_base in Pa - Psi atmospheric = 101325 Pa at sea level
    Jv in m**3/s
    length in m

    """
    # topology
    #n=40

    g = markov.linear(n)
    g = compute_flux(g,n=n, psi_e=psi_e, psi_base=psi_base, Jv=Jv, k0=k0, length=length)
    #scene = plot(g, prop_cmap='J_out', has_radius=True)
    scene= None
    return g, scene


def test_tree(g = None, n=30, psi_e=300000, psi_base=101325, Jv=1e-10, k0=0.5, length=1e-4, prop_cmap='radius'):
    """ Test flux and water potential computation on a linear root. 
    
    Units :
    psi_e & psi_out in Pa
    Jv in m**3/s
    length in m

    """
    # topology
    #n=40
    if g is None :
        g = markov.markov_binary_tree(nb_vertices=n)
    g = compute_flux(g,n=n, psi_e=psi_e, psi_base=psi_base, Jv=Jv, k0=k0, length=length)
    scene = plot(g, prop_cmap=prop_cmap, has_radius=True)
    return g, scene

def root_visitor(g, v, turtle):
    angles = [90,45]+[30]*5
    n = g.node(v)
    radius = n.radius*1.e4
    order = n.order
    length = n.length*1.e4

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


def plot(g =None, length=1.e-4, has_radius=False, r_base=1., r_tip=0.25, visitor=root_visitor, prop_cmap='radius'):
    """
    Exemple:

        >>> from openalea.plantgl.all import *
        >>> s = plot()
        >>> shapes = dict( (x.getId(), x.geometry) for x in s)
        >>> Viewer.display(s)
    """
    if g is None:
        g = random_binary_tree()

    # compute length
    for v in g:
        n = g.node(v)
        n.length = length 
    
    if not has_radius:
        discont_radius(g,r_base=r_base, r_tip=r_tip)

    turtle = turt.PglTurtle()
    turtle.down(180)
    scene = turt.TurtleFrame(g, visitor=root_visitor, turtle=turtle, gc=False)

    # Compute color from radius
    my_colormap(g,prop_cmap)

    shapes = dict( (sh.getId(),sh) for sh in scene)

    colors = g.property('color')
    for vid in colors:
        shapes[vid].appearance = pgl.Material(colors[vid])
    scene = pgl.Scene(shapes.values())
    return scene

def my_colormap(g, property_name, cmap='jet',lognorm=True):
    prop = g.property(property_name)
    keys = prop.keys()
    values = np.array(prop.values())
    m, M = int(values.min()), int(values.max())
    print m, M
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
    

