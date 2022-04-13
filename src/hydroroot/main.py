

from warnings import warn
import numpy as np
from hydroroot.length import fit_law
from hydroroot import radius, flux, conductance 
from hydroroot.generator import markov, measured_root # 21-12-14: FB __init__.py in src not doing job


def hydroroot_mtg(
    primary_length=0.15,
    delta=2.e-3,
    beta=0.25,
    order_max=5,
    segment_length=1e-4,
    nude_length=0.02,
    seed=2,
    ref_radius=1e-4,
    order_decrease_factor=0.7,
    length_data=None,
    n=None,
    **kwds
):
    """
    Simulate and generate a root system.
    
    Parameters
    ----------
        - primary_length (Float) the primary root length
        - delta (Float) the inter-branching length
        - beta (Float) branching variability, random variability around delta (0.25 gives 25% of variability)
        - order_max (int) root order maximum
        - segment_length (Float) the MTG segment length
        - nude_length (Float) distance to the tip without any laterals
        - seed (int) the seed for the random generator
        - ref_radius (Float) the primary root radius
        - order_decrease_factor (Float) the radius decrease factor applied when increasing order
        - length_data: (pandas dataframe) the length law 
        - n: (int) maximum number of vertices

    Returns
    -------
        - g
        - surface
        - volume

    Example
    =======

    """
    xl, yl = length_data
    length_law = fit_law(xl, yl, scale=segment_length)

    # compute the architecture
    nb_nude_vertices = int(nude_length / segment_length)
    branching_delay = int(delta / segment_length)

    if n is None:
        nb_vertices = int(primary_length / segment_length)
    else:
        nb_vertices = n
        warn("Use primary_length instead")

    g = markov.markov_binary_tree(
        nb_vertices=nb_vertices,
        branching_variability=beta,
        branching_delay=branching_delay,
        length_law=length_law,
        nude_tip_length=nb_nude_vertices,
        order_max=order_max,
        seed=seed)

    # compute radius property on MTG
    g = radius.ordered_radius(g, ref_radius=ref_radius, order_decrease_factor=order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    return g, surface, volume

def hydroroot_flow(
    g,
    segment_length=1.e-4,
    k0=300,
    Jv=0.1,
    psi_e=0.4,
    psi_base=0.1,
    axial_conductivity_data=None,
    radial_conductivity_data=None,
):
    """
    Flux and equivalent conductance calculation

    Parameters
    ----------
        - g (MTG)
        - segment_length (Float) the MTG segment length
        - k0 (Float) not used
        - Jv (Float) flux in (microL/s)
        - psi_e (Float) external hydrostatic potential (MPa)
        - psi_base (Float) hydrostatic potential at the root base (MPa)
        - axial_conductivity_data (list of Float) K vs distance to tip
        - radial_conductivity_data (list of Float) k vs distance to tip

    Returns
    -------
	- g (MTG): the MTG with 

    """
    xa, ya = axial_conductivity_data
    # commented line below, BUG correction, the global flux was diverging when decreasing segment_length
    # ya was not supposed to be multiplied by (segment_length / 1.e-4)
    # ya = list(np.array(ya) * (segment_length / 1.e-4))
    axial_conductivity_law = fit_law(xa, ya)

    xr, yr = radial_conductivity_data
    radial_conductivity_law = fit_law(xr, yr)

    # Compute K using axial conductance data
    g = conductance.fit_property_from_spline(g, axial_conductivity_law, 'position', 'K_exp')
    g = conductance.compute_K(g) # Fabrice 2020-01-17: calculation of K in dimension [L^3 P^(-1) T^(-1)]
    # Compute the flux

    g = conductance.fit_property_from_spline(g, radial_conductivity_law, 'position', 'k0')
    g = conductance.compute_k(g, k0='k0')

    # TODO: return Keq base and Jv
    g = flux.flux(g, Jv, psi_e, psi_base, invert_model=True)

    Keqs = g.property('Keq')
    v_base = next(g.component_roots_at_scale_iter(g.root, scale=1))

    Keq = Keqs[v_base]
    Jv_global = Keq * (psi_e - psi_base)

    return g, Keq, Jv_global


def hydroroot(
    primary_length=0.15,
    delta=2.e-3,
    beta=0.25,
    order_max=5,
    segment_length=1e-4,
    nude_length=0.02,
    seed=2,
    ref_radius=1e-4,
    order_decrease_factor=0.7,
    k0=300,
    Jv=0.1,
    psi_e=0.4,
    psi_base=0.1,
    length_data=None,
    axial_conductivity_data=None,
    radial_conductivity_data=None,
    n=None
):
    """Simulate a root system and compute global conductance and flux.

    see hydroroot_mtg and hydroroot_flow

    """
    g, surface, volume = hydroroot_mtg(primary_length=primary_length,
                                       delta=delta,
                                       beta=beta,
                                       order_max=order_max,
                                       segment_length=segment_length,
                                       nude_length=nude_length,
                                       seed=seed,
                                       ref_radius=ref_radius,
                                       order_decrease_factor=order_decrease_factor,
                                       length_data=length_data,
                                       n=n,
                                       )
    xa, ya = axial_conductivity_data
    # commented line below, BUG correction, the global flux was diverging when decreasing segment_length
    # ya was not supposed to be multiplied by (segment_length / 1.e-4)
    # ya = list(np.array(ya) * (segment_length / 1.e-4))
    axial_conductivity_law = fit_law(xa, ya)

    xr, yr = radial_conductivity_data
    radial_conductivity_law = fit_law(xr, yr)

    # Compute K using axial conductance data
    g = conductance.fit_property_from_spline(g, axial_conductivity_law, 'position', 'K_exp')
    g = conductance.compute_K(g)  # Fabrice 2020-01-17: calculation of K in dimension [L^3 P^(-1) T^(-1)]
    # Compute the flux

    g = conductance.fit_property_from_spline(g, radial_conductivity_law, 'position', 'k0')
    g = conductance.compute_k(g, k0='k0')

    # TODO: return Keq base and Jv
    g = flux.flux(g, Jv, psi_e, psi_base, invert_model=True)

    Keqs = g.property('Keq')
    v_base = next(g.component_roots_at_scale_iter(g.root, scale=1))

    Keq = Keqs[v_base]
    Jv_global = Keq * (psi_e - psi_base)

    return g, surface, volume, Keq, Jv_global


def hydroroot_from_data(
    primary_length=0.15,
    delta=2.e-3,
    beta=0.25,
    order_max=5,
    segment_length=1e-4,
    nude_length = 0.02,
    seed = 2,
    ref_radius = 1e-4,
    order_decrease_factor = 0.7,
    k0 = 300,
    Jv = 0.1,
    psi_e = 0.4,
    psi_base = 0.1,
    length_data=None,
    axial_conductivity_data=None,
    radial_conductivity_data=None,
    primary_length_data=None,
    lateral_length_data=None,
    ):
    """ Reconstruct a root system and compute global conductance and flux.

    Parameters
    ==========
    - primary_length_data (list Float) data with the lateral positions from the base
    - lateral_length_data (list Float) data with the lateral lengths

    see hydroroot_mtg and hydroroot_flow for the other parameters

    """
    xl, yl = length_data
    length_law = fit_law(xl, yl, scale=segment_length)

    xa, ya = axial_conductivity_data
    # commented line below, BUG correction, the global flux was diverging when decreasing segment_length
    # ya was not supposed to be multiplied by (segment_length / 1.e-4)
    # ya = list(np.array(ya) * (segment_length / 1.e-4))
    axial_conductivity_law = fit_law(xa, ya)

    xr, yr = radial_conductivity_data
    radial_conductivity_law = fit_law(xr, yr)

    # compute the architecture
    nb_nude_vertices = int(nude_length / segment_length)
    branching_delay = int(delta / segment_length)

    g = measured_root.mtg_builder(
        primary_length,
        primary_length_data,
        lateral_length_data,
        segment_length=segment_length,
        branching_variability=beta,
        branching_delay=branching_delay,
        length_law=length_law,
        nude_tip_length=nb_nude_vertices,
        order_max=order_max,
        seed=seed)


    # compute radius property on MTG
    g = radius.ordered_radius(g, ref_radius=ref_radius, order_decrease_factor=order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # Compute K using axial conductance data
    g = conductance.fit_property_from_spline(g, axial_conductivity_law, 'position', 'K_exp')
    g = conductance.compute_K(g)  # Fabrice 2020-01-17: calculation of K in dimension [L^3 P^(-1) T^(-1)]

    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    # Compute the flux

    g = conductance.fit_property_from_spline(g, radial_conductivity_law, 'position', 'k0')
    g = conductance.compute_k(g, k0='k0')

    # TODO: return Keq base and Jv
    g = flux.flux(g, Jv, psi_e, psi_base, invert_model=True)

    Keqs = g.property('Keq')
    v_base = next(g.component_roots_at_scale_iter(g.root, scale=1))

    Keq = Keqs[v_base]
    Jv_global = Keq * (psi_e - psi_base)

    return g, surface, volume, Keq, Jv_global
