

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
    """Simulate and generate a root system.

    :param primary_length: Float (Default value = 0.15)
    :param delta: Float (Default value = 2.e-3)
    :param beta: Float (Default value = 0.25)
    :param order_max: int (Default value = 5)
    :param segment_length: Float (Default value = 1e-4)
    :param nude_length: Float (Default value = 0.02)
    :param seed: int (Default value = 2)
    :param ref_radius: Float (Default value = 1e-4)
    :param order_decrease_factor: Float (Default value = 0.7)
    :param length_data: pandas dataframe (Default value = None)
    :param n: int (Default value = None)

    :returns: - g
        - surface
        - volume
    
    Example

    """
    # F. Bauget 2022-08-12: added if-else to be able to use the function without length data, useful for usage demo
    if length_data:
        xl, yl = length_data
        length_law = fit_law(xl, yl, scale=segment_length)
    else:
        length_law = None

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
    """Flux and equivalent conductance calculation

    :param g: MTG
    :param segment_length: Float (Default value = 1.e-4)
    :param k0: Float (Default value = 300)
    :param Jv: Float (Default value = 0.1)
    :param psi_e: Float (Default value = 0.4)
    :param psi_base: Float (Default value = 0.1)
    :param axial_conductivity_data: list of Float (Default value = None)
    :param radial_conductivity_data: list of Float (Default value = None)
    :returns: - g (MTG): the MTG with the following properties filled: K (axial conductance), k (radial donductivity),
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
    see :func:`~main.hydroroot_mtg` and  :func:`~main.hydroroot_flow`

    :param primary_length:  (Default value = 0.15)
    :param delta:  (Default value = 2.e-3)
    :param beta:  (Default value = 0.25)
    :param order_max:  (Default value = 5)
    :param segment_length:  (Default value = 1e-4)
    :param nude_length:  (Default value = 0.02)
    :param seed:  (Default value = 2)
    :param ref_radius:  (Default value = 1e-4)
    :param order_decrease_factor:  (Default value = 0.7)
    :param k0:  (Default value = 300)
    :param Jv:  (Default value = 0.1)
    :param psi_e:  (Default value = 0.4)
    :param psi_base:  (Default value = 0.1)
    :param length_data:  (Default value = None)
    :param axial_conductivity_data:  (Default value = None)
    :param radial_conductivity_data:  (Default value = None)
    :param n:  (Default value = None)

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
    """Reconstruct a root system and compute global conductance and flux.

    :param primary_length_data: list Float (Default value = None)
    :param lateral_length_data: list Float (Default value = None)
    :param return: 
    :param g: MTG
    :param surface: float
    :param volume: float
    :param Keq: float
    :param Jv_global: float
    :param see: func
    :param primary_length:  (Default value = 0.15)
    :param delta:  (Default value = 2.e-3)
    :param beta:  (Default value = 0.25)
    :param order_max:  (Default value = 5)
    :param segment_length:  (Default value = 1e-4)
    :param nude_length:  (Default value = 0.02)
    :param seed:  (Default value = 2)
    :param ref_radius:  (Default value = 1e-4)
    :param order_decrease_factor:  (Default value = 0.7)
    :param k0:  (Default value = 300)
    :param Jv:  (Default value = 0.1)
    :param psi_e:  (Default value = 0.4)
    :param psi_base:  (Default value = 0.1)
    :param length_data:  (Default value = None)
    :param axial_conductivity_data:  (Default value = None)
    :param radial_conductivity_data:  (Default value = None)

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
    """wrapper function: build a MTG with properties like radius and vertex length set.
    The MTG is either generated or built from a data.
    The radius and vertex length properties are set.
    The following properties are computed: length, position, mylength, surface, volume, total length, primary root length, nb of intercepts

    :param primary_length: primary root length for generated mtg (Default value = 0.13)
    :param seed: seed for generated mtg (Default value = None)
    :param delta: branching delay for generated mtg (Default value = 2.0e-3)
    :param nude_length: length from tip without lateral for generated mtg (Default value = 2.0e-2)
    :param df: pandas DataFrame with the architecture data to be reconstructed (Default value = None)
    :param segment_length: float (Default value = 1.0e-4)
    :param length_data: string (Default value = None)
    :param branching_variability: float (Default value = 0.25)
    :param order_max: float (Default value = 4.0)
    :param order_decrease_factor: float (Default value = 0.7)
    :param ref_radius: float (Default value = 7.0e-5)
    :param intercepts: list
    :param Flag_radius: boolean (Default value = False)
    :returns: - g: MTG with the different properties set or computed (see comments above),
        - primary_length: primary root length (output for generated mtg)
        - total_length: total root length
        - surface: total root surface
        - intercepts: nb of intercepts at a given distance from base
        - _seed: the seed used in the generator

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
