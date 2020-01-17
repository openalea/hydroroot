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
Flux computation.
"""


from openalea.mtg import traversal
from openalea.mtg.traversal import *
import numpy as np

class Flux(object):   # edit this to also allow for flux computation instead just redistribution
    """Compute the water potential and fluxes at each vertex of the MTG.

    """

    def __init__(self, g,
                 Jv, psi_e, psi_base,
                 invert_model=False,
                 k=None, K=None, CONSTANT=1.,
                 cut_and_flow=False):
        """ Flux computes water potential and fluxes at each vertex of the MTG `g`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `k` (dict) - lateral conductance
            - `K` (dict) - axial conductance
            - `Jv` (float) - water flux at the root base in microL/s
            - `psi_e` - hydric potential outside the roots (pressure chamber) in MPa
                if None, then consider that the value has been defined on each vertex.
            - `psi_base` - hydric potential at the root base (e.g. atmospheric pressure for decapited plant) in MPa
            - `invert_model` - when false, distribute output flux within the root ; when true, compute the output flux
                for the given root and conditions.
            - `cut_and_flow` (bool) - Use specific model to compute conductance at tips with cut & flow.
        :Example:

            flux = Flux(g, ...)
        """
        self.CONSTANT = CONSTANT  # used for sensitivity analysis
        self.g = g
        self.k = k if k else g.property('k')
        self.K = K if K else g.property('K')
        self.Jv = Jv
        self.psi_e = psi_e if psi_e else g.property('psi_e')
        self.psi_base = psi_base
        self.length = g.property('length')
        self.invert_model = invert_model

        self.HAS_SOIL = psi_e is None
        self.CUT_AND_FLOW = cut_and_flow

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
        # Fabrice 2020-01-17: commented lines below because a MTG property must not changed see conductance.compute_K
        #                   and here K is a shadow copy of g.property('K') then g.property('K') is modified below
        # for vid in K:
        #     K[vid] /= length[vid]
        #     K[vid] *= CONSTANT

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
            kids = g.children(v)
            if (not self.CUT_AND_FLOW) or kids:
                r = 1./(k[v] + sum(Keq[cid] for cid in kids))
                R = 1./K[v]
                Keq[v] = 1./(r+R)
            else:
                assert self.CUT_AND_FLOW and (len(kids)==0)
                Keq[v] = K[v]
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

                if not self.HAS_SOIL:
                    psi_in[v] = (K[v] * psi_out[v] + psi_e * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                    j[v] = (psi_e - psi_in[v]) * k[v]
                else:
                    psi_in[v] = (K[v] * psi_out[v] + psi_e[v] * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                    j[v] = (psi_e[v] - psi_in[v]) * k[v]

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
                if not self.HAS_SOIL:
                    psi_in[v] = (K[v] * psi_out[v] + psi_e * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                else:
                    psi_in[v] = (K[v] * psi_out[v] + psi_e[v] * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)

            #print 'exiting Psi computation'

            print 'entering Jv computation'
            for v in traversal.post_order2(g, v_base):
            # compute water flux according to the psis from root tips to root base
            # modif Fabrice 2020-01-17: if CUT_AND_FLOW then the input conductance at the tip is the axial one
            #     kin = K[v] if self.CUT_AND_FLOW and len(g.children(v)) == 0 else k[v]
                kin = K[v] if self.CUT_AND_FLOW and g.nb_children(v) == 0 else k[v]
                if not self.HAS_SOIL:
                    # j[v] = (psi_e - psi_in[v]) * k[v]
                    j[v] = (psi_e - psi_in[v]) * kin
                else:
                    # j[v] = (psi_e[v] - psi_in[v]) * k[v]
                    j[v] = (psi_e[v] - psi_in[v]) * kin
            # End modif Fabrice 2020-01-17
                children = g.children_iter(v)
                # if children is None: #Fabrice 2020-01-17: wrong syntax never None even if there are no children use len instead
                if len(g.children(v)) == 0:
                    J_out[v] = j[v]
                else:  # TODO CHECK THIS !!!
                    influx = j[v] + sum( J_out[cid] for cid in children )
                    J_out[v] = influx #(psi_in[v]-psi_out[v])*K[v]
                #print J_out

            if not self.HAS_SOIL:
                Jv_global = Keq[v_base] * (psi_e - psi_base)
            else:
                Jv_global = Keq[v_base] * (psi_e[v_base] - psi_base)

            print "Local Computation Water Flux Jvl = ", J_out[v_base]
            print "Global Computation Water Flux Jvg = ", Jv_global



class RadialShuntFlux(Flux):
    """ Compute the water potential and fluxes at each vertex of the MTG.

    On each vertex, the topology of the radial resistance network
    has one direct shortcut to the parent. The shortcut has two new parameters a and b.
    """

    def __init__(self, a=1., b=0., **kwds):
        """ Flux computes water potential and fluxes at each vertex of the MTG `g`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `k` (dict) - lateral conductance
            - `K` (dict) - axial conductance
            - `Jv` (float) - water flux at the root base in microL/s
            - `psi_e` - hydric potential outside the roots (pressure chamber) in MPa
                if None, then consider that the value has been defined on each vertex.
            - `psi_base` - hydric potential at the root base (e.g. atmospheric pressure for decapited plant) in MPa
            - `invert_model` - when false, distribute output flux within the root ; when true, compute the output flux for the given root and conditions
            - `a` - relative factor to the main radial path conductivity.
            - `b` - relative factor to the shortcut path conductivity.

        :Algorithm:

        :Example:

            flux = RadialShuntFlux(g, ..., a=0.8, b=0.2)
        """
        Flux.__init__(self, **kwds)
        self.a = a if a else g.property('a')
        self.b = b if b else g.property('b')


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

        # NEW addition
        a = self.a
        b = self.b

        # Select the base of the root
        v_base = g.component_roots_at_scale_iter(g.root, scale=g.max_scale()).next()

        # Add properties
        g.add_property('Keq')
        g.add_property('psi_in')
        g.add_property('psi_out')
        g.add_property('j')
        g.add_property('J_out')

        # NEW properties
        g.add_property('alpha')
        g.add_property('beta')

        # Convert axial conductivities to axial conductances
        # Fabrice 2020-01-17: commented lines below because a MTG property must not changed see conductance.compute_K
        #                   and here K is a shadow copy of g.property('K') then g.property('K') is modified below
        # for vid in K:
        #     K[vid] /= length[vid]
        #     K[vid] *= CONSTANT

        # Equivalent conductance computation
        Keq = g.property('Keq')
        alpha = g.property('alpha')
        beta = g.property('beta')

        #print 'entering Keq computation'

        # TODO
        # NEW EQUATION g, children Keq, k, K, a, b
        # equation
        for v in traversal.post_order2(g, v_base):
            # compute
            r, ra, rb, R = 0., 0., 0., 0.
            # compute Keq[v]
            #r = 1./(k[v] + sum(Keq[cid] for cid in g.children_iter(v)))
            #R = 1./K[v]
            #Keq[v] = 1./(r+R)
            # compute alpha and beta
            Keqc = sum(Keq[cid] for cid in g.children_iter(v))
            alpha[v] = a*b*k[v] + (a+b)*K[v]+b*Keqc
            alpha[v] /= a*b*k[v] + K[v]*(1+a+b)

            beta[v] = K[v] + b*Keqc
            beta[v] /= a*b*k[v] + K[v]*(1+a+b)

        #print 'exiting Keq computation'

        # Water flux and water potential computation
        psi_out = g.property('psi_out')
        psi_in = g.property('psi_in')
        j = g.property('j')
        J_out = g.property('J_out')



        # TODO Write the two systems
        if not(invert_model) :
            # distribute a given output into the root system

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

                if not self.HAS_SOIL:
                    # TODO
                    #psi_in[v] = (K[v] * psi_out[v] + psi_e * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                    #j[v] = (psi_e - psi_in[v]) * k[v]
                    pass
                else:
                    # TODO
                    #psi_in[v] = (K[v] * psi_out[v] + psi_e[v] * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                    #j[v] = (psi_e[v] - psi_in[v]) * k[v]
                    pass

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

                # TODO
                if not self.HAS_SOIL:
                    psi_in[v] = (K[v] * psi_out[v] + psi_e * (-a * k[v] * beta[v] + Keq_children)) / (a * k[v] + K[v] + Keq_children - a * k[v] * alpha[v])
                    #psi_in[v] = (K[v] * psi_out[v] + psi_e * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)
                else:
                    psi_in[v] = (K[v] * psi_out[v] + psi_e[v] * (-a * k[v] * beta[v] + Keq_children)) / (a * k[v] + K[v] + Keq_children - a * k[v] * alpha[v])
                    # psi_in[v] = (K[v] * psi_out[v] + psi_e[v] * (k[v] + Keq_children)) / (k[v] + K[v] + Keq_children)

            #print 'exiting Psi computation'

            #print 'entering Jv computation'
            for v in traversal.post_order2(g, v_base):
            # compute water flux according to the psis from root tips to root base
                if not self.HAS_SOIL:
                    j[v] = (psi_e - psi_in[v]) * k[v] * alpha[v]
                else:
                    j[v] = (psi_e[v] - psi_in[v]) * k[v] * alpha[v]
                children = g.children_iter(v)
                if children is None:
                    J_out[v] = j[v]
                else:  # TODO CHECK THIS !!!
                    influx = j[v] + sum( J_out[cid] for cid in children )
                    J_out[v] = influx #(psi_in[v]-psi_out[v])*K[v]

            if not self.HAS_SOIL:
                Jv_global = Keq[v_base] * (psi_e - psi_base)
            else:
                Jv_global = Keq[v_base] * (psi_e[v_base] - psi_base)


def flux(g, Jv=0.1, psi_e=0.4, psi_base=0.101325,
         invert_model=False, k=None, K=None, CONSTANT=1.,
         shunt=False, a=1., b=0.,
         cut_and_flow=False):
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
            - `shunt` (bool) : use the RadialShunt Flux (True) or the classical one (False)
            - `a` : relative factor to the main radial path conductivity.
            - `b` : relative factor to the shortcut path conductivity.

        :Example::

            my_flux = flux(g)
    """
    if not shunt:
        f = Flux(g, Jv, psi_e, psi_base, invert_model, k=k, K=K, CONSTANT=CONSTANT, cut_and_flow=cut_and_flow)
    else:
        f = RadialShuntFlux(a, b, g=g, Jv=Jv, psi_e=psi_e, psi_base=psi_base, invert_model=invert_model, k=k, K=K, CONSTANT=CONSTANT)

    f.run()

    return f.g

def segments_at_length(g, l, root=1, dl=1e-4):
    """Returns all the segments intercepted at a given length.

    Parameters
    ==========
        - g: MTG
        - l: length

    Returns
    =======
        - number of segment
    """
    length = {}

    if 'mylength' in g.property_names():
        length = g.property('mylength')
    else:
        for v in traversal.pre_order2(g, root):
            pid = g.parent(v)
            length[v] = length[pid] + dl if pid else dl
        g.properties()['mylength'] = length

    vids = []
    for v in g:
        pid = g.parent(v)
        if pid and (length[pid] <= l <= length[v]):
            vids.append(v)
    return vids


def cut(g, cut_length, threshold=1e-4):
    # Added Fabrice 2020-01-17: segment_length in parameters list
    """Cut the architecture at a given length `cut_length`.

        :Parameters:
            - `g` (MTG) - the root architecture
            - `cut_length` (float, m) - length at which the architecture is cut from collar.
            - 'segment_length' (float, mm) - length of the vertices

        :Returns:
            - `g`(MTG) - the architecture after the cut process. This is a copy.

        :Example::

            g_cut = cut(g, 0.09) # Cut g at 9cm. Remove the 2 last cm of a root architecture of 11 cm (primary length).
    """
    # vids = segments_at_length(g, cut_length)
    vids = segments_at_length(g, cut_length, dl = threshold)

    g_cut = g.copy()
    for v in vids:
        # g_cut.remove_tree(v)
        # the for loop below is a copy of openalea.mtg.Tree.remove_tree but use post_order2 instead of post_order
        #    to avoid "RuntimeError: maximum recursion depth exceeded"
        for vtx_id in post_order2(g, v):
            g_cut.remove_vertex(vtx_id)

    return g_cut


def ramification_length_law(g, root=1, dl=1e-4):
    """Returns the length of the ramified axes along the main axis.

    X is the distance to tip along the main axis
    Y is the length of the ramification
    """
    length = {}

    if 'back_length' in g.property_names():
        length = g.property('back_length')
    else:
        for v in traversal.post_order2(g, root):
            sids = g.Sons(v, EdgeType='<')
            length[v] = length[sids[0]] + dl if sids else dl
        g.properties()['back_length'] = length


    axis1 = g.Axis(root)

    ramif_length = []
    total_length = length[root]
    for v in axis1[:-1]:
        if g.nb_children(v) > 1:
            ramification_ids = g.Sons(v, EdgeType='+')
            for rid in ramification_ids:
                ramif_length.append((length[v], length[rid]))
    ramif_length = list(reversed(ramif_length))

    X, Y = zip(*ramif_length)
    X = np.array(X)

    X /= total_length
    Y = np.array(Y)

    return X, Y
