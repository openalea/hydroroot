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


class Flux(object):   # edit this to also allow for flux computation instead just redistribution
    """ Compute the water potential and fluxes at each vertex of the MTG.

    """

    def __init__(self, g, Jv, psi_e, psi_base, invert_model=False, k=None, K=None, CONSTANT=1.):
        """ Flux computes water potential and fluxes at each vertex of the MTG `g`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `k` (dict) - lateral conductance
            - `K` (dict) - axial conductance
            - `Jv` (float) - water flux at the root base in microL/s
            - `psi_e` - hydric potential outside the roots (pressure chamber) in MPa
            - `psi_base` - hydric potential at the root base (e.g. atmospheric pressure for decapited plant) in MPa
            - `invert_model` - when false, distribute output flux within the root ; when true, compute the output flux for the given root and conditions

        :Example:

            flux = Flux(g, ...)
        """
        self.CONSTANT = CONSTANT  # used for sensitivity analysis
        self.g = g
        self.k = k if k else g.property('k')
        self.K = K if K else g.property('K')
        self.Jv = Jv
        self.psi_e = psi_e
        self.psi_base = psi_base
        self.length = g.property('length')
        self.invert_model = invert_model

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
        
        g = self.g; k = self.k; K = self.K ; CONSTANT = self.CONSTANT
        Jv = self.Jv; psi_e = self.psi_e; psi_base = self.psi_base
        length = self.length; invert_model = self.invert_model

        # Select the base of the root
        v_base = g.component_roots_at_scale_iter(g.root, scale=g.max_scale()).next()

        # Add properties
        g.add_property('Keq')
        g.add_property('psi_in')
        g.add_property('psi_out')
        g.add_property('j')
        g.add_property('J_out')

        # Convert axial conductivities to axial conductances
        for vid in K:
            K[vid] /= length[vid]
            K[vid] *= CONSTANT

        # Apply scaling k and K values
        #for vid in k:
        #    k[vid] *= CONSTANT
        #for vid in K:
        #    K[vid] *= CONSTANT
        #Jv *= CONSTANT

        # Equivalent conductance computation
        Keq = g.property('Keq')
        #print 'entering Keq computation'
        for v in traversal.post_order2(g, v_base):
            r = 1./(k[v] + sum(Keq[cid] for cid in g.children_iter(v)))
            R = 1./K[v]
            Keq[v] = 1./(r+R)
        #print 'exiting Keq computation'

        # Water flux and water potential computation
        psi_out = g.property('psi_out')
        psi_in = g.property('psi_in')
        j = g.property('j')
        J_out = g.property('J_out')

        if not(invert_model) : # distribute a given output into the root system

            #print 'entering Jv distribution'

            for v in traversal.pre_order2(g, v_base):
            #compute psi according to Millman theorem, then compute radial flux
                parent = g.parent(v)
                brothers = g.children_iter(parent)
                children = g.children_iter(v)

                Keq_brothers = sum( Keq[cid] for cid in brothers)
                Keq_children = sum( Keq[cid] for cid in children)

                if parent is None:
                    assert v == v_base
                    psi_out[v] = psi_base
                    J_out[v] = Jv
                else:
                    psi_out[v] = psi_in[parent]
                    J_out[v] = (J_out[parent] - j[parent]) * ( Keq[v] / Keq_brothers )

                psi_in[v] = (K[v] * psi_out[v] + psi_e * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                j[v] = (psi_e - psi_in[v]) * k[v]

            #print 'exiting Jv distribution'

        else :  # compute the water output for the given root system and conditions

            #print 'entering Psi computation'
            for v in traversal.pre_order2(g, v_base):
            #compute psi according to Millman theorem from root base to root tips
                parent = g.parent(v)
                children = g.children_iter(v)
                if parent is None:
                    assert v == v_base
                    psi_out[v] = psi_base
                else:
                    psi_out[v] = psi_in[parent]
                Keq_children = sum( Keq[cid] for cid in children )
                psi_in[v] = (K[v] * psi_out[v] + psi_e * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
            #print 'exiting Psi computation'

            #print 'entering Jv computation'
            for v in traversal.post_order2(g, v_base):
            # compute water flux according to the psis from root tips to root base
                j[v] = (psi_e - psi_in[v]) * k[v]
                children = g.children_iter(v)
                if children is None:
                    J_out[v] = j[v]
                else:  # TODO CHECK THIS !!!
                    influx = j[v] + sum( J_out[cid] for cid in children )
                    J_out[v] = influx #(psi_in[v]-psi_out[v])*K[v]

            Jv_global = Keq[v_base]*(psi_e-psi_base)

            #print 'exiting Jv computation'

            #print "Keq base = ", Keq[v_base]
            #print "psi_e, psi_base = ", psi_e, psi_base
            #print "Local Computation Water Flux Jvl = ", J_out[v_base]
            #print "Global Computation Water Flux Jvg = ", Jv_global


#        for v in traversal.pre_order2(g, v_base):
#        # this version of the algorithm causes numerical issues - do not use
#            parent = g.parent(v)
#            if parent is None:
#                assert v == v_base
#                psi_out[v] = psi_base
#                #print 'psi_out',v,psi_out[v]
#                J_out[v] = Jv
#                #print 'j_out',v,J_out[v]
#            else:
#                psi_out[v] = psi_in[parent]
#                #print 'psi_out',v,psi_out[v]
#                J_out[v] = (J_out[parent] - j[parent]) * ( Keq[v] / (sum( Keq[cid] for cid in g.children(parent))))
#                #print 'j_out',v,J_out[v]
#            psi_in[v] = (J_out[v] / K[v]) + psi_out[v]
#            #print 'psi_in',v,psi_in[v]
#            j[v] = (psi_e-psi_in[v]) * k[v]
#            #print 'j',v,j[v]

#            if J_out[v] < j[v]:
#                print 'Vertex %d (Jout=%.4f, j=%.4f, psi_in=%.4f)'%(v,J_out[v]/CONSTANT, j[v]/CONSTANT,psi_in[v]/CONSTANT)

        # UnNormalize k and K values
#X         for vid in k:
#X             k[vid] /= CONSTANT
#X         for vid in K:
#X             K[vid] /= CONSTANT
#X         for vid in Keq:
#X             Keq[vid] /= CONSTANT
#X         for vid in j:
#X             j[vid] /= CONSTANT
#X         for vid in J_out:
#X             J_out[vid] /= CONSTANT

        #Jv *= CONSTANT


def flux(g, Jv=0.1, psi_e=0.4, psi_base=0.101325, invert_model=False, k=None, K=None, CONSTANT=1.):
    """ flux computes water potential and fluxes at each vertex of the MTG `g`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `Jv` (float) - water flux at the root base in microL/s
            - `psi_e` - hydric potential outside the roots (pressure chamber) in MPa
            - `psi_base` - hydric potential at the root base (e.g. atmospheric pressure for decapited plant) in MPa
            - `invert_model` - when false, distribute output flux within the root ; when true, compute the output flux for the given root and conditions


        :Optional Parameters:
            - `k` (dict) - lateral conductance
            - `K` (dict) - axial conductance

        :Example::

            my_flux = flux(g)
    """
    f = Flux(g, Jv, psi_e, psi_base, invert_model, k=k, K=K, CONSTANT=CONSTANT)
    f.run()

    return f.g
