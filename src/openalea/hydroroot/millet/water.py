from __future__ import absolute_import

# CPL: Not used. Remove it.
# from warnings import warn
# import numpy as np
from openalea.hydroroot.length import fit_law
from openalea.hydroroot import radius, markov, flux, conductance
from openalea.hydroroot.millet import architecture,law
from math import log


def hydrolic(g,
             segment_length=1e-4,
             ref_radius=1e-4,
             order_decrease_factor=0.7,
             k0=300,
             Jv=0.1,
             psi_e=0.4,
             psi_base=0.1,
             length_data=None,
             axial_conductivity_data=None,
             radial_conductivity_data=None,
             ):
    """ Simulate a root system and compute global conductance and flux.

    Parameters
    ==========

    Returns
    =======
        - surface
        - volume
        - Keq
        - Jv

    Example
    =======

    """

    HAS_SOIL = psi_e is None

    xl, yl = length_data
    length_law = fit_law(xl, yl, scale=segment_length)

    # xa, ya = axial_conductivity_data
    # ya = list(np.array(ya) * (segment_length / 1e-4))
    # axial_conductivity_law = fit_law(xa, ya)

    xr, yr = radial_conductivity_data
    radial_conductivity_law = fit_law(xr, yr)

    # compute the architecture
    # nb_nude_vertices = int(nude_length / segment_length)
    # branching_delay = int(delta / segment_length)

    # compute radius property on MTG
    # TODO: Add a different set of parameters for crown and seminal roots.
    # g = radius.ordered_radius(g, ref_radius=ref_radius, order_decrease_factor=order_decrease_factor)

    # CPL: import modules outside the body of functions
    # import math

    # CPL: Define these laws outside the function.
    seminal_law = lambda x: 70.42 * log(x) + 282.45
    crown_law = lambda x: 59.15 * log(x) + 416.57
    lateral_law = lambda x: 29.95 * log(x) + 286.32
    g = law.compute_diameter(g,seminal_RootDiameter_law=seminal_law,crown_RootDiameter_law=crown_law, lr_RootDiameter_law=lateral_law)
    g = law.radius_from_computed_diameters(g)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # Compute K using axial conductance data
    # g = conductance.fit_property_from_spline(g, axial_conductivity_law, 'position', 'K')
    g = conductance.compute_K_from_laws(g)

    # Compute surface and volume
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    # Compute the flux

    g = conductance.fit_property_from_spline(g, radial_conductivity_law, 'position', 'k0')
    g = conductance.compute_k(g, k0='k0')

    # TODO: return Keq base and Jv
    g = flux.flux(g, Jv, psi_e, psi_base, invert_model=True)

    Keqs = g.property('Keq')
    v_base = g.component_roots_at_scale_iter(g.root, scale=1).next()
    Keq = Keqs[v_base]

    if HAS_SOIL:
        psi_e_base = g.property('psi_e')[v_base]
    else:
        psi_e_base = psi_e
    Jv_global = Keq * (psi_e_base - psi_base)

    return g, surface, volume, Keq, Jv_global
