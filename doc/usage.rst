=====
Usage
=====

Many examples are given in the notebook *example\boursiac2022.ipynb* that reproduces most of the figures and tables of Boursiac et al. [boursiac2022]_.

Run calculation on a given architecture
---------------------------------------
The following example shows how to run a simulation to calculate the out going sap flux and the equivalent conductance of a given architecture. Hydroroot is able to read two file format for architecture: RSML (http://rootsystemml.github.io/) and a simple tabulated text file that only details root lengths and branching positions. The later format is enough to describe root in hydroponic condition.

The text format architecture file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=======================  ========================  =====
distance_from_base_(mm)  lateral_root_length_(mm)  order
=======================  ========================  =====
0.89                     90.81             	       1
3.02                     63.98             	       1
102.94                     0.0             	       1
2.14                     23.72             	       1-1
90.81                     0.0             	       1-1
2.48                     5.15             	       1-2
63.98                     0.0             	       1-2
=======================  ========================  =====

This is a tab separated text file with 3 columns:

1. the distance from base of the branching laterals in mm
2. the lateral root length in mm
3. a string of one or more number indicating the parent root

In the example above, the root has two lateral of 1st order and on each of them one lateral of 2d order. The order of 1 indicates that the laterals are on the primary root. The last line with order 1 with 0.0 in the second column indicates the primary root tip.
The line with order 1-1 indicates that this a second order lateral on the first lateral positioned at 2.14 mm from the branching on the primary root. And so on.

A 4th column with the averaged diameter of the root may be given, that may be used to build the MTG representing the architecture.

Running the calculation
~~~~~~~~~~~~~~~~~~~~~~~

The corresponding notebook is *example/example_archi_from_file.ipynb*

If the package HydroRoot is not installed, the following examples can be run by cloning the sources from git and then sourcing the src directory in Ipython console for instance like this:

.. code-block:: python

    import sys
    sys.path.extend(['../src']) # if run from the example folder for instance

assuming the dependencies installed.

The following lines present a small example of calculation of the sap flux from an Arabidopsis de-topped root plunged in a hydroponic solution at a hydrostatic pressure of 0.4 Mpa when its  base is at the atmospheric pressure.

.. code-block:: python

    from hydroroot.display import plot
    from hydroroot.read_file import read_archi_data
    from hydroroot.main import hydroroot_flow, root_builder

Reading the file architecture as a DataFrame

.. code-block:: python

    df = read_archi_data('data/plant-01.txt')

Building the MTG from the file, and  return the primary root length, the total length and the surface. The seed refer to the seed of the root generator when the MTG is not built from a file but is generated.

.. code-block:: python

    g, primary_length, total_length, surface, seed = root_builder(df=df, segment_length=1.0e-4)

Some conductance data versus distance to tip

.. code-block:: python

    k_radial_data=([0, 0.2],[30.0,30.0])
    K_axial_data=([0, 0.2],[3.0e-7,4.0e-4])

Flux and equivalent conductance calculation, for a root in an external hydroponic medium at 0.4 MPa, its base at 0.1 MPa,
and with the conductances set above.

.. code-block:: python

    g, keq, jv = hydroroot_flow(g, psi_e = 0.4, psi_base = 0.1, axial_conductivity_data = K_axial_data, radial_conductivity_data = k_radial_data)

.. code-block:: python

    print('equivalent root conductance (microL/s/MPa): ',keq, 'sap flux (microL/s): ', jv)

Displaying the water uptake along the architecture using the Plantgl viewer (https://github.com/openalea/plantgl).

.. code-block:: python

    %gui qt
    plot(g, prop_cmap='j') # j is the radial influx in ul/s

You may change the property to display the hydrostatic pressure inside the xylem vessels for instance

.. code-block:: python

    plot(g, prop_cmap='psi_in') # P in MPa

You may change the radial conductivity and see the impact on the water uptake

.. code-block:: python

    k_radial_data=([0, 0.2],[300.0,300.0])
    g, keq, jv = hydroroot_flow(g, psi_e = 0.4, psi_base = 0.1, axial_conductivity_data = K_axial_data, radial_conductivity_data = k_radial_data)
    print('sap flux (microL/s): ', jv)
    plot(g, prop_cmap='j')

Or the axial conductance

.. code-block:: python

    k_radial_data=([0, 0.2],[30.0,30.0])
    K_axial_data=([0, 0.2],[3.0e-7,1.0e-4])
    g, keq, jv = hydroroot_flow(g, psi_e = 0.4, psi_base = 0.1, axial_conductivity_data = K_axial_data, radial_conductivity_data = k_radial_data)
    print('sap flux (microL/s): ', jv)
    plot(g, prop_cmap='j')

Importing architecture from RSML
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is a small example to illustrate how to use the RSML format (http://rootsystemml.github.io/). The architecture is the
arabidopsis-simple example http://rootsystemml.github.io/images/examples/arabidopsis-simple.rsml.

.. code-block:: python

    import rsml
    from hydroroot import radius
    from hydroroot.main import hydroroot_flow
    from hydroroot.display import plot
    from hydroroot.hydro_io import import_rsml_to_discrete_mtg, export_mtg_to_rsml

We first read the RSML file and convert it into a *continuous* MTG. This is a MTG where each root (primary and lateral)
is represented by one vertex. The geometry of each root is then stored in g_c.property(‘geometry’).

.. code-block:: python

    g_c = rsml.rsml2mtg('data/arabidopsis-simple.rsml')

To be used in HydroRoot the MTG has to be converted to a *discrete* form of MTG, i.e. each vertex represent a representative
elementary volume of a given length for example :math:`10^{-4}` m. In HydroRoot the lengths are in meter, therefore we
must retrieve the resolution and unit of the RSML file,

.. code-block:: python

    resolution = g_c.graph_properties()['metadata']['resolution'] # pixel to unit
    unit = g_c.graph_properties()['metadata']['unit']
    print(unit)

In this example, the resolution of the RSML file is 0.01 and the unit is cm. The length unit in HydroRoot is the meter.
Therefore, to pass from pixels (the raw data in the RSML file) to the meter we must multiply *g_c.graph_properties()['metadata']['resolution']*
by 0.01.

.. code-block:: python

    resolution = resolution * 0.01 # pixel to unit to m

We build the discrete MTG.

.. code-block:: python

    g = import_rsml_to_discrete_mtg(g_c, segment_length = 1.0e-4, resolution = resolution)

We calculate some properties needed to simulate a sap flux from the root.

.. code-block:: python


    g = radius.ordered_radius(g, 7.0e-5, 0.7) # root radii
    g = radius.compute_relative_position(g) # Compute the position of each segment relative to the axis bearing it

Some conductance data versus distance to tip

.. code-block:: python

    k_radial_data=([0, 0.2],[30.0,30.0])
    K_axial_data=([0, 0.2],[3.0e-7,4.0e-4])

Flux and equivalent conductance calculation, for a root in an external hydroponic medium at 0.4 MPa, its base at 0.1 MPa,
and with the conductances set above.

.. code-block:: python

    g, keq, jv = hydroroot_flow(g, psi_e = 0.4, psi_base = 0.1, axial_conductivity_data = K_axial_data, radial_conductivity_data = k_radial_data)

Display the local water uptake heatmap in 3D

.. code-block:: python

    %gui qt
    plot(g, prop_cmap = 'j')

We may also export the MTG to RSML

.. code-block:: python

    export_mtg_to_rsml(g, "test.rsml", segment_length = 1.0e-4)

The resolution of the RSML data is 1.0e-4 and the unit is meter.
At this stage (2022-08-22) only the root length and the branching
position are used to simulate architecture in hydroponic solution. The
exact position in 3D is not stored in the discrete MTG form and so not
exported to RMSL.

Run calculation on a generated architecture
-------------------------------------------

The corresponding notebook is *example/example_generated_archi.ipynb*

If the examples are run using the source, add the source directory to the system path

.. code-block:: python

    import sys;
    sys.path.extend(['../src'])

.. code-block:: python

    import pandas 
    from hydroroot.main import root_builder, hydroroot_flow
    from hydroroot.display import plot

The Hydroroot generator of architecture is described in Boursiac et al. [boursiac2022]_.
It uses length distribution law for laterals, specific to a given species, to generate realistic architecture. Here we use the length laws determined for Col0 arabidopsis.

.. code-block:: python

    length_data = [] # length law used to generate arabidopsis realistic architecture
    for filename in ['data/length_LR_order1_160615.csv','data/length_LR_order2_160909.csv']:
        df = pandas.read_csv(filename, sep = ';', header = 1, names = ('LR_length_mm', 'relative_distance_to_tip'))
        df.sort_values(by = 'relative_distance_to_tip', inplace = True)
        length_data.append(df)

We generate the MTG with some specific parameters: 
 + primary_length:length of the primary root 
 + delta: the average distance between lateral branching 
 + branching_variability: the variability of the branching distance around delta 
 + nude_length: distance from the tip without any laterals 
 + order_max: the maximum order of laterals

And return the primary root length, the total length and the surface.  Seed may be used as seed to generate the same architecture.

.. code-block:: python

    g, primary_length, total_length, surface, seed = root_builder(primary_length = 0.13, delta = 2.0e-3, nude_length = 2.0e-2, segment_length = 1.0e-4,
                                                      length_data = length_data, branching_variability = 0.25, order_max = 4.0, order_decrease_factor = 0.7,
                                                      ref_radius = 7.0e-5)



Some conductance data versus distance to tip

.. code-block:: python

    k_radial_data=([0, 0.2],[30.0,30.0])
    K_axial_data=([0, 0.2],[3.0e-7,4.0e-4])

Flux and equivalent conductance calculation, for a root in an external
hydroponic medium at 0.4 MPa, its base at 0.1 MPa, and with the
conductances set above.

.. code-block:: python

    g, keq, jv = hydroroot_flow(g, psi_e = 0.4, psi_base = 0.1, axial_conductivity_data = K_axial_data, radial_conductivity_data = k_radial_data)

.. code-block:: python

    print(keq,jv)


.. parsed-literal::

    0.007146429180199128 0.002143928754059739


Display the local water uptake heatmap in 3D

.. code-block:: python

    %gui qt
    plot(g, prop_cmap='j') # j is the radial flux in ul/s

Model parameters
----------------

The main model parameters are grouped in the python class parameters, see :class:`hydroroot.init_parameter.Parameters`.
The parameters may be passed to the class by reading a yaml file, see :meth:`hydroroot.init_parameter.Parameters.read_file`.

There are two solvers in HydroRoot project. The first, used for the paper Boursiac et al. 2022 [boursiac2022]_, is a purely water transport model. The second is a solute and water transport model. Therefore, the *solute* category in the yaml file has meaning only for the second solver.

The structure of the yaml file is the following (see examples at https://github.com/openalea/hydroroot)

| **archi**
|	**read_architecture:** Boulean
|		True read an architecture file, False generate an architecture
|	**input_dir**: String
|		the folder with the architecture file, relative path to the script
|	**input_file**: list of string
|		list of architecture file names, eg. [file1.txt] or [file1.txt, file2.txt, file3.txt] wildcar may be used
|	**seed**: int or list of int
|		the seed used to generate architecture
|	**length_file**: list of string
|		name of the files containing the length law, relative path
|		file format: "LR_length_mm" ; "relative_distance_to_tip"
|		laws used to generate lateral roots of the 1st order (1_order_law), and lateral roots of order above 1 (2_order_law)
|	**primary_length**: float or list of float
|		primary root length in m used when generating architecture
|		unit: m
|	**branching_delay**: float or list of float
|		distance between branching
|		unit: m
|	**branching_variability**: float
|		maximum random variation around the branching_delay value
|		between [0 ; 1]
|	**order_max**: int
|		maximum order of laterals possible
|	**segment_length**: float
|		MTG vertices length
|		unit: m
|	**nude_length**: float or list of float
|		part of roots without any lateral root, distance from tip
|		unit: m
|	**ref_radius**: float
|		reference radius of the primary root
|		unit: m
|	**order_decrease_factor**: float
|		radius decrease factor applied when increasing order
|		radius of lateral of order n: :math:`r = \beta^n R_{ref}`
|		with :math:`r = \beta` order_decrease_factor and :math:`R_{ref}` ref_radius
| **hydro**
|	**k0**: float
|		radial conductivity
|		unit: :math:`\mu L.s^{-1}.MPa^{-1}.m^{-2}`
|	**axial_conductance_data**: 2 list of float
|		axial conductance versus distance to tip, K(x)
|		unit: :math:`\mu L.m.s^{-1}.MPa^{-1}`
| **solute**
|  **J_s**: float
|   	active pumping rate
|   	unit: mol/(m2.s)
|  **P_s**: float
|		permeability coefficient
|		unit: m/s
|  **Cse**: float
|   	concentration of permeating solutes
|       unit: :math:`mol.m^{-3} \text{or}\ mM`
|  **Ce**: float
|   	concentration of non-permeating solutes
|       unit: :math:`mol.m^{-3} \text{or}\ mM`
|  **sigma**: float
|   	reflection coefficient
|   	dimensionless
| **experimental**
|	**Jv**:  float
|		flux at the root base
|		unit: :math:`\mu L.s^{-1}`
|	**psi_e**:  float
|		hydrostatic pressure outside the root (pressure chamber)
|		unit: :math:`MPa`
|	**psi_base**:  float
|		hydrostatic pressure at the root base (e.g. atmospheric pressure for decapitated plant)
|		unit: :math:`MPa`
| **output**:
|	**intercepts**: float or list of float
|		distance from the base at which the number of intercepts are calculated
|		unit: m
|	**radfold**: float or list of float
|		factor to explore a k0 range
|	**axfold**: float or list of float
|		factor to explore a axial conductance range
|	**run_nb**: int
|		number of run with the same set of parameters

Few parameters may be set to lists allowing to run successive simulations.
For list of number there are two syntax: [x1, ..., xn] or range(start, end, step).
For example, range(0.02, 0.09, 0.02) or [0.02, 0.04, 0.06, 0.08] will give the same results.
The parameter will take successively the values 0.02, 0.04, 0.06 and 0.08.
The parameter *run_nb*  would be useful with read_architecture = False and no given seed to generate different architectures.

**Note:** Parameter is just a python class. It can not be used directly with Hydroroot functions, intermediary script should be used.
We will give you some examples using scripts that be found at https://github.com/openalea/hydroroot in example.

Run simple calculation using the Parameters class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The corresponding notebook is *example/example_parameter_class.ipynb*

.. code-block:: python

    import sys; print('Python %s on %s' % (sys.version, sys.platform))
    sys.path.extend(['../src'])


.. parsed-literal::

    Python 3.8.12 | packaged by conda-forge | (default, Jan 30 2022, 23:42:07) 
    [GCC 9.4.0] on linux


.. code-block:: python

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

.. code-block:: python

    parameter = Parameters()
    parameter.read_file('parameters_palnt_01.yml')

Read the architecture file and build the MTG

.. code-block:: python

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

.. code-block:: python

    g, Keq, Jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                psi_e = parameter.exp['psi_e'],
                                psi_base = parameter.exp['psi_base'],
                                axial_conductivity_data = parameter.hydro['axial_conductance_data'],
                                radial_conductivity_data = parameter.hydro['k0'])

.. code-block:: python

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
    

.. code-block:: python

    plot(g, prop_cmap='j') # j is the radial flux in ul/s

Example of solute and water transport simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example

.. code-block:: python

    import sys; print('Python %s on %s' % (sys.version, sys.platform))
    sys.path.extend(['../src'])


.. parsed-literal::

    Python 3.8.12 | packaged by conda-forge | (default, Jan 30 2022, 23:42:07) 
    [GCC 9.4.0] on linux


.. code-block:: python

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

.. code-block:: python

    parameter = Parameters()
    parameter.read_file('parameters_Ctr-3P2.yml')

In the code the concentration are in :math:`mol.\mu L^{-1}`

.. code-block:: python

    Cse = parameter.solute['Cse'] * 1e-9 # mol/m3 -> mol/microL, external permeating solute concentration
    Ce = parameter.solute['Ce'] * 1e-9 # mol/m3 -> mol/microL, external non-permeating solute concentration

Read the architecture file and build the MTG

.. code-block:: python

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

.. code-block:: python

    g = set_conductances(g, axial_pr = parameter.hydro['axial_conductance_data'], k0_pr = parameter.hydro['k0']) 
    g = flux.flux(g, psi_e = parameter.exp['psi_e'], psi_base = parameter.exp['psi_base'])  # initialization
    g = init_some_MTG_properties(g, tau = parameter.solute['J_s'], Cini = parameter.solute['Cse'])

Perform the calculation, this a Newtown-Raphson loop on a matrix system,
then there is a convergence loop.

.. code-block:: python

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

.. code-block:: python

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

.. code-block:: python

    plot(g, prop_cmap='C') # C is the radial flux in mol/microL
