from math import pi
from collections import defaultdict

from openalea.mtg import *
from openalea.mtg import algo

def poiseuille(radius, length, viscosity=1e-3):
    """
    Compute a conductance of a xylem element based on their radius and length.
    
    Parameters
    ==========
    radius : float (m)
        radius of a xylem tube

    length: float (m) 
        length of a xylem element

    viscosity : float (Pa.s)
        dynamic viscosity of the liquid
    
    The poiseuille formula is:
        :math:` conductance = \frac{\pi r^4}{8 \mu L }` 
        with :math:`r` the radius of a pipe, 
        :math:`\mu` the viscosity of the liquid,
        :math:`L` the length of the pipe.
        
    .. seealso:: http://en.wikipedia.org/wiki/Poiseuille
    """
    return pi*(radius**4) / ( 8 * viscosity * length)


def compute_k(g, k0 = 0.1, length=1.e-4):
    """ Set radial conductances (k) in a MTG at a given value. """
    radius = g.property('radius')
    k = dict( (vid,radius[vid]*2*pi*length*k0) for vid in g.vertices(scale=g.max_scale()))
    return k


def compute_K(g, length=1.e-4, nb_xylem=5, radius_scale = 1/10.):
    """ 
    Set axial conductances (K) in a MTG according to Poiseuille law. 

    The conductance depends on the radius of each xylem pipe, the number of xylem pipes,
    and on the length of a root segment.

    radius_scale allows to compute the radius of a xylem pipe from the radius of a root segment.
    """

    radius = g.property('radius_xylem')
    if not radius:
        full_radius = g.property('radius')
        radius = dict( (vid,r*radius_scale) for vid,r in full_radius.iteritems())
    nb_xylem = g.property('nb_xylem')
    if not nb_xylem:
        nb_xylem = defaultdict(lambda : 5)
    K = dict((vid, nb_xylem[vid]*poiseuille(radius[vid], length)) for vid in g.vertices(scale=g.max_scale()))
    return K


