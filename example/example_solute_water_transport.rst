Example

.. code:: ipython3

    import sys; print('Python %s on %s' % (sys.version, sys.platform))
    sys.path.extend(['../src'])


.. parsed-literal::

    Python 3.8.12 | packaged by conda-forge | (default, Jan 30 2022, 23:42:07) 
    [GCC 9.4.0] on linux


.. code:: ipython3

    import math
    from hydroroot import flux
    from hydroroot.main import root_builder
    from hydroroot.init_parameter import Parameters
    from hydroroot.display import plot
    from hydroroot.read_file import read_archi_data
    from hydroroot.conductance import set_conductances
    from hydroroot.water_solute_transport import pressure_calculation_no_non_permeating_solutes, init_some_MTG_properties
    
    # for the PlantGL viewer used in hydroroot.display.plot
    %gui qt 

Read the yaml file and set the Parameters variables, assuming that the
code is run from the example folder

.. code:: ipython3

    parameter = Parameters()
    parameter.read_file('parameters_Ctr-3P2.yml')

In the code the concentration are in :math:`mol.\mu L^{-1}`

.. code:: ipython3

    Cse = parameter.solute['Cse'] * 1e-9 # mol/m3 -> mol/microL, external permeating solute concentration
    Ce = parameter.solute['Ce'] * 1e-9 # mol/m3 -> mol/microL, external non-permeating solute concentration

Read the architecture file and build the MTG

.. code:: ipython3

    fname = parameter.archi['input_dir'] + parameter.archi['input_file'][0]
    df = read_archi_data(fname)
    g, primary_length, total_length, surface, seed = root_builder( primary_length = parameter.archi['primary_length'],
                                                                    delta = parameter.archi['branching_delay'],
                                                                    nude_length = parameter.archi['nude_length'], 
                                                                    df = df,
                                                                    segment_length = parameter.archi['segment_length'],
                                                                    length_data = parameter.archi['length_data'],
                                                                    order_max = parameter.archi['order_max'],
                                                                    order_decrease_factor = parameter.archi['order_decrease_factor'],
                                                                    ref_radius = parameter.archi['ref_radius'])

Set the conductance in the MTG (in previous examples that was done in
hydroroot_flow), set some other properties in *init_some_MTG_properties*
and perform some initialization. Note that here *parameter.hydro[‘k0’]*
is a float.

.. code:: ipython3

    g = set_conductances(g, axial_pr = parameter.hydro['axial_conductance_data'], k0_pr = parameter.hydro['k0']) 
    g = flux.flux(g, psi_e = parameter.exp['psi_e'], psi_base = parameter.exp['psi_base'])  # initialization
    g = init_some_MTG_properties(g, tau = parameter.solute['J_s'], Cini = parameter.solute['Cse'])

Perform the calculation, this a Newtown-Raphson loop on a matrix system,
then there is a convergence loop.

.. code:: ipython3

    eps = 1.0e-9 # global: stop criterion for the Newton-Raphson loop in Jv_P_calculation and Jv_cnf_calculation
    nb_v = g.nb_vertices()
    Fdx = 1.0
    Fdx_old = 1.
    while Fdx > eps:
        g, dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g, sigma = parameter.solute['Sigma'], 
                                                                               tau = parameter.solute['J_s'], 
                                                                               Ce = Ce,
                                                                               Ps = parameter.solute['P_s'], 
                                                                               Cse = Cse, 
                                                                               Pe = parameter.exp['psi_e'], 
                                                                               Pbase = parameter.exp['psi_base'])
        Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
        if abs(Fdx - Fdx_old) < eps: break
        Fdx_old = Fdx
    Jv = g.property('J_out')[1]

.. code:: ipython3

    result=f"""
    primary length (m): {primary_length}
    surface (m2): {surface}
    total length (m): {total_length}
    flux (microL/s): {Jv}
    """
    print(result)


.. parsed-literal::

    
    primary length (m): 0.434
    surface (m2): 0.005643500494241343
    total length (m): 3.979
    flux (microL/s): 0.025700314390202567
    


Display the concentration in the architecture

.. code:: ipython3

    plot(g, prop_cmap='C') # C is the radial flux in mol/microL

