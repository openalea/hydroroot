import random
import numpy as np

from matplotlib import pyplot, mpl
from pylab import cm, colorbar
from pylab import plot as pylab_plot
from matplotlib.colors import Normalize, LogNorm

from openalea.mtg import turtle as turt
import openalea.plantgl.all as pgl
from openalea.deploy.shared_data import shared_data

import hydroroot
from hydroroot.flux import *
from hydroroot import radius, conductance, markov, display
from hydroroot.read_file import readCSVFile

def compute_flux(g, n=300, psi_e=400000., psi_base=101325., Jv=1e-10, k0=0.3e-12, length=1e-4):
    k0 = float(k0)
    g = radius.discont_radius(g, r_base=1.e-4, r_tip=5.e-5)
    g = radius.compute_length(g,length)
    g = radius.compute_relative_position(g)
    surf = radius.compute_surface(g)
    print('surf',surf)
    g = conductance.compute_k(g, k0)
    g = compute_conductance(g)

    k = g.property('k')
    K = g.property('K')
    assert all(v>0 for v in list(K.values())),K

    g = flux(g, Jv, psi_e, psi_base)

    J_out = g.property('J_out')
    #assert all(v>0 for v in J_out.values()), J_out.values()
    return g

def compute_conductance(g, fn='conductance_data.csv'):
    fn = shared_data(hydroroot, fn, share_path='share')
    data = readCSVFile(fn)
    g = conductance.fit_property_from_csv(g, data, 'position', 'K', k=1)
    return g

def test_linear(n=1500, psi_e=400000., psi_base=101325., Jv=1e-10, k0=0.3e-12, length=1e-4):
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


def test_tree(g = None, n=1500, lr=0.9,  psi_e=400000, psi_base=101325, Jv=1e-10, k0=0.3e-12, length=1e-4, prop_cmap='radius'):
    """ Test flux and water potential computation on a linear root.

    Units :
    psi_e & psi_out in Pa
    Jv in m**3/s
    length in m

    """
    # topology
    #n=40
    if g is None :
        g = markov.markov_binary_tree(nb_vertices=n, branching_stability=lr, seed=2)
    g = compute_flux(g,n=n, psi_e=psi_e, psi_base=psi_base, Jv=Jv, k0=k0, length=length)
    scene = display.plot(g, prop_cmap=prop_cmap, has_radius=True)
    return g, scene


