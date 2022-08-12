

from warnings import warn
# import numpy as np

from openalea.mtg import traversal

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
	- g (MTG): the MTG with the following properties filled: K (axial conductance), k (radial donductivity),
	            j (radial flux), J_out (axial flux), psi_in and psi_out (hydrostatic pressure into the root at the
	            input and output of a MTG node
	- Keq (float): the equivalent conductance of the whole root
	- Jv_global (float): the outgoing flux at the root base

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

    :parameters:
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

def root_builder(primary_length = 0.13, seed = None, delta = 2.0e-3, nude_length = 2.0e-2, df = None, segment_length = 1.0e-4,
                  length_data = None, branching_variability = 0.25, order_max = 4.0, order_decrease_factor = 0.7,
                  ref_radius = 7.0e-5, Flag_radius = False):
    """
    wrapper function: build a MTG with properties like radius and vertex length set.

    The MTG is either generated or built from a data.
    The radius and vertex length properties are set.
    The following properties are computed: length, position, mylength, surface, volume, total length,
        primary root length, nb of intercepts

    :Parameters:
        primary_length: primary root length for generated mtg
        seed:  seed for generated mtg, if None randomly generated
        delta: branching delay  for generated mtg
        nude_length: length from tip without lateral for generated mtg
        df: pandas DataFrame with the architecture data to be reconstructed
        segment_length: float (1.0e-4) - MTG segment length
        length_data: string - length laws file names
        branching_variability: float (0.25) - random variability of delta
        order_max: float (4.0) - maximum lateral order
        order_decrease_factor: float (0.7) - radius decrease factor between order
        ref_radius: float (7.0e-5) - the primary root radius
        intercepts: list (None) - list of distance from base
        Flag_radius: boolean (False) - False: radius calculated in radius.ordered_radius; True: radius from the architecture file
    :Returns:
        g: MTG with the different properties set or computed (see comments above),
        primary_length: primary root length (output for generated mtg)
        total_length: total root length
        surface: total root surface
        intercepts: nb of intercepts at a given distance from base
        _seed: the seed used in the generator
    """
    if df is not None:
        g = measured_root.mtg_from_aqua_data(df, segment_length)
        _seed = None
    else:
        # if no seed just create one
        if seed is None:
            _seed = markov.my_seed()
        else:
            _seed = seed

        g = markov.generate_g(_seed, length_data,
                       branching_variability, delta,
                       nude_length, primary_length, segment_length,
                       order_max)

    # compute radius property on MTG
    # F. Bauget 2022-05-17 : added if
    if (not Flag_radius) or ('radius' not in g.property_names()):
        g = radius.ordered_radius(g, ref_radius, order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # Calculation of the distance from base of each vertex, used for cut and flow
    _mylength = {}
    for v in traversal.pre_order2(g, 1):
        pid = g.parent(v)
        _mylength[v] = _mylength[pid] + segment_length if pid else segment_length
    g.properties()['mylength'] = _mylength

    # total_length is the total length of the RSA (sum of the length of all the segments)
    total_length = g.nb_vertices(scale = 1) * segment_length
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    if df is not None:
        v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))
        primary_length = g.property('position')[v_base]

    return g, primary_length, total_length, surface, _seed
