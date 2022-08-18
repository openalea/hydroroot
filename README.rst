OpenAlea.HydroRoot
==================

Software
~~~~~~~~~~~~

Authors
-------
  * Christophe Pradal
  * Yann Boursiac
  * Fabrice Bauget
  * Mikael Lucas
  * Christophe Godin
  * Christophe Maurel

Institutes
----------
CIRAD / INRA / inria / CNRS

Status
------
Python package

License
-------
CecILL-C

**URL** : http://openalea.rtfd.io

About
~~~~~~

Description
-----------

OpenAlea.HydroRoot is a hydraulic root architecture modelling and a root architecture system generator package.


Content
-------

The OpenAlea.HydroRoot package contains a root architecture simulation model coupled with an hydraulic solver or
a solute-water transportsolver.

Example
-------

Heat map representaion of the incoming local radial flows on an arabidopsis root.

.. image:: example/data/fig-6E.png
    :width: 600

Installation
------------
::

    conda create -n hydroroot3 -c conda-forge -c openalea3 openalea.hydroroot
    conda activate hydroroot3

Requirements
------------
  - numpy
  - openalea.deploy
  - openalea.mtg
  - pandas >0.17
  - path
  - python >=3.6
  - pyyaml
  - rsml
  - scipy

Usage
-----

The following example builds a root and display it in the PlantGL viewer.

.. code-block:: python

    from hydroroot.main import root_builder
    from hydroroot.display import plot
    %gui qt
    g, primary_length, total_length, surface, seed = root_builder(order_max=1)
    plot(g)


The following example builds a root, run the hydraulic solver and display the eat map representaion of the incoming
local radial flows on an arabidopsis root in the PlantGL viewer.

K is the axial conductance in $10^{-9} m^4.s^{-1}.MPa^{-1}$, k is radial conductivity in $10^{-9} m.s^{-1}.MPa^{-1}$

.. code-block:: python

    from hydroroot.main import hydroroot
    from hydroroot.display import plot
    %gui qt
    K = ([0,0.2],[0.0,1.0e-2])
    k = ([0.0,0.2],[300.0,300.0])
    g, surface, volume, Keq, Jv_global = hydroroot(axial_conductivity_data = K, radial_conductivity_data=k, order_max = 1)
    plot(g, prop_cmap = 'j')

See also the jupyter notebook boursiac2022.ipynb in example/Bourisac2022 for examples. This notebook is aimed to run different simulations to
generate figures and tables of Bourisac et al. 2022 (see below), illustrating the HydroRoot capabilities.

Citation
--------
Yann Boursiac, Christophe Pradal, Fabrice Bauget, Mikaël Lucas, Stathis Delivorias, Christophe Godin, Christophe Maurel,
Phenotyping and modeling of root hydraulic architecture reveal critical determinants of axial water transport,
Plant Physiology, 2022;, kiac281, https://doi.org/10.1093/plphys/kiac281

Documentation
~~~~~~~~~~~~~
https://hydroroot.rtfd.io

