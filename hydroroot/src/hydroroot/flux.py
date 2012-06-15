# -*- python -*-
#
#       HydroRoot
#
#       Copyright 2012 CNRS - INRIA - CIRAD - INRA  
#
#       File author(s): Mikael Lucas <mikael.lucas.at.supagro.inra.fr>
#                       Christophe Pradal <christophe.pradal.at.cirad.fr>
#                       Christophe Maurel 
#                       Christophe Godin
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
################################################################################

"""

"""

from openalea.mtg import traversal

class Flux(object):
    """ Compute the water potential and fluxes at each vertex of the MTG.

    """

    def __init__(self, g, Jv, psi_e, psi_base, k=None, K=None):
        """ Flux computes water potential and fluxes at each vertex of the MTG `g`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `k` (dict) - lateral conductance
            - `K` (dict) - axial conductance
            - `Jv` (float) - water flux at the root base in m**3/s
            - `psi_e` - hydric potential outside the roots (pressure chamber) in Pa
            - `psi_base` - hydric potential at the root base (e.g. atmospheric pressure for decapited plant) in Pa

        :Example:

            flux = Flux(g, ...)
        """
        self.g = g
        self.k = k if k else g.property('k')
        self.K = K if K else g.property('K')
        self.Jv = Jv
        self.psi_e = psi_e
        self.psi_base = psi_base

    def run(self):
        """ Compute the water potential and fluxes of each segments

        For each vertex of the root, compute :
            - the water potential (:math:`\psi^{\text{out}}`) at the base;
            - the water flux (`J`) at the base;
            - the lateral water flux (`j`) entering the segment.

        :Algorithm:
            The algorithm has two stages:
                - First, on each segment, an equivalent conductance is computed in post_order (children before parent).
                - Finally, the water flux and potential are computed in pre order (parent then children).
        """
        
        g = self.g; k = self.k; K = self.K
        Jv = self.Jv; psi_e = self.psi_e; psi_base = self.psi_base

        # Select the base of the root
        v_base = g.component_roots_at_scale(g.root, scale=g.max_scale()).next()

        # Add properties
        g.add_property('Keq')
        g.add_property('psi_in')
        g.add_property('psi_out')
        g.add_property('j')
        g.add_property('J_out')

        # Conductance computation
        Keq = g.property('Keq')
        for v in traversal.post_order2(g, v_base):
            r = 1./(k[v] + sum(Keq[cid] for cid in g.children(v))) 
            R = 1./K[v]
            Keq[v] = 1./(r+R)

        # Water flux and water potential computation
        psi_out = g.property('psi_out')
        psi_in = g.property('psi_in')
        j = g.property('j')
        J_out = g.property('J_out')

        for v in traversal.pre_order2(g, v_base):
            parent = g.parent(v)
            if parent is None:
                assert v == v_base
                psi_out[v] = psi_base
                #print 'psi_out',v,psi_out[v]
                J_out[v] = Jv
                #print 'j_out',v,J_out[v]
            else:
                psi_out[v] = psi_in[parent]
                #print 'psi_out',v,psi_out[v]
                J_out[v] = ((J_out[parent] - j[parent]) * Keq[v]) / (sum( Keq[cid] for cid in g.children(parent)))
                #print 'j_out',v,J_out[v]

            psi_in[v] = (J_out[v] / K[v]) + psi_out[v]
            #print 'psi_in',v,psi_in[v]
            j[v] = (psi_e-psi_in[v]) * k[v]
            #print 'j',v,j[v]


def flux(g, Jv=1e-10, psi_e=300000., psi_base=101325., k=None, K=None):
    """ flux computes water potential and fluxes at each vertex of the MTG `g`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `Jv` (float) - water flux at the root base in m**3/s
            - `psi_e` - hydric potential outside the roots (pressure chamber) in Pa
            - `psi_base` - hydric potential at the root base (e.g. atmospheric pressure for decapited plant) in Pa


        :Optional Parameters:
            - `k` (dict) - lateral conductance
            - `K` (dict) - axial conductance

        :Example::

            my_flux = flux(g)
    """    
    f = Flux(g, Jv, psi_e, psi_base, k=k, K=K)
    f.run()
    return f.g
