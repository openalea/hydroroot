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
from hydroroot import radius, conductance, markov, display



def surface(g, length=1e-4):
    surf = 0
    radius = g.property('radius')
    for vid in g.vertices():
        if radius.has_key(vid):
            surf += 2 * pi * radius[vid] * length
    return surf

def compute_flux(g, n=300, psi_e=300000., psi_base=101325., Jv=1e-10, k0=0.0005, length=1e-4):
    k0 = float(k0)
    radius.discont_radius(g, r_base=1.e-4, r_tip=5.e-5)
    surf = surface(g)
    print 'surf',surf
    k = conductance.compute_k(g, k0,length)
    K = conductance.compute_K(g,length=length)

    assert all(v>0 for v in K.values()),K

    g = flux(g, k, K, Jv, psi_e, psi_base)

    J_out = g.property('radius')
    #assert all(v>0 for v in J_out.values()), J_out.values()
    return g

def test_linear(n=300, psi_e=300000., psi_base=101325., Jv=1e-10, k0=0.005, length=1e-4):
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
    scene = display.plot(g, prop_cmap='J_out', has_radius=True)
    return g, scene


def test_tree(g = None, n=30, psi_e=300000, psi_base=101325, Jv=1e-10, k0=0.005, length=1e-4, prop_cmap='radius'):
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
    scene = display.plot(g, prop_cmap=prop_cmap, has_radius=True)
    return g, scene

