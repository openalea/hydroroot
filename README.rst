OpenAlea.HydroRoot 
==================

Software
~~~~~~~~~~~~

Authors
-------
  * Christophe Pradal
  * Yann Boursiac
  * Mikael Lucas
  * Fabrice Bauget
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

The OpenAlea.HydroRoot package contains a root architecture simulation model coupled with an hydraulic solver. 

Example
-------

Heat map representaion of the incoming local radial flows on an arabidopsis root. 
![Alt Text](example/data/fig-6E.png)

Installation
------------

Conda Installation
++++++++++++++++++
::

    conda create -n hydroroot -c conda-forge -c openalea3 openalea.hydroroot


Requirements 
++++++++++++

    * openalea.mtg
    * openalea.plantgl
    * openalea.visualea
    * RSML
    * pandas > 0.17
    * numpy
    * scipy

Usage
+++++
::

    see the jupyter notebook figures_tables.ipynb for examples. This notebook is aimed to run different simulations to
    generate figures and tables illustrating the HydroRoot capabilities.

Documentation
~~~~~~~~~~~~~
https://hydroroot.rtfd.io

Citations
~~~~~~~~~

If you use Hydroroot for your research, please cite:

1. Yann Boursiac, Christophe Pradal, Fabrice Bauget, Mikaël Lucas, Stathis Delivorias, Christophe Godin, Christophe Maurel, Phenotyping and modeling of root hydraulic architecture reveal critical determinants of axial water transport, Plant Physiology, Volume 190, Issue 2, October 2022, Pages 1289–1306, https://doi.org/10.1093/plphys/kiac281

2. Fabrice Bauget, Virginia Protto, Christophe Pradal, Yann Boursiac, Christophe Maurel, A root functional-structural model allows to assess effects of water deficit on water and solute transport parameters, Journal of Experimental Botany, 2022;, erac471, https://doi.org/10.1093/jxb/erac471

