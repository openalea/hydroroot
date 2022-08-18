Example

.. code:: ipython3

    import sys; print('Python %s on %s' % (sys.version, sys.platform))
    sys.path.extend(['../src'])


.. parsed-literal::

    Python 3.8.12 | packaged by conda-forge | (default, Jan 30 2022, 23:42:07) 
    [GCC 9.4.0] on linux


.. code:: ipython3

    import pandas as pd
    from hydroroot import radius
    from hydroroot.main import hydroroot_flow, root_builder
    from hydroroot.init_parameter import Parameters
    from hydroroot.generator.measured_root import mtg_from_aqua_data
    from hydroroot.display import plot
    from hydroroot.read_file import read_archi_data
    
    # for the PlantGL viewer used in hydroroot.display.plot
    %gui qt 

Read the yaml file and set the Parameters variables, assuming that the
code is run from the example folder

.. code:: ipython3

    parameter = Parameters()
    parameter.read_file('parameters_palnt_01.yml')

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

Calculation of the equivalent conductance and the sap flux

.. code:: ipython3

    g, Keq, Jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                psi_e = parameter.exp['psi_e'],
                                psi_base = parameter.exp['psi_base'],
                                axial_conductivity_data = parameter.hydro['axial_conductance_data'],
                                radial_conductivity_data = parameter.hydro['k0'])

.. code:: ipython3

    result=f"""
    primary length (m): {primary_length}
    surface (m2): {surface}
    total length (m): {total_length}
    flux (microL/s): {Jv}
    """
    print(result)


.. parsed-literal::

    
    primary length (m): 0.10300000000000001
    surface (m2): 0.0004625701757655344
    total length (m): 1.6260000000000001
    flux (microL/s): 0.0028789143185531108
    


.. code:: ipython3

    plot(g, prop_cmap='j') # j is the radial flux in ul/s

