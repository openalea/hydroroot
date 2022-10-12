============
Installation
============

First install miniconda following the instructions given here https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

User installation
---------------------

HydroRoot may by installed simply on a conda environments:

::

    conda create -n hydroroot -c conda-forge -c openalea3 openalea.hydroroot
    conda activate hydroroot

That creates a conda environment called *hydroroot*, install in it *openalea.hydroroot* with all the dependencies and
activate the environment. Then just open an Ipython session and enjoy.

If you want notebook support, run for example:

::

    conda install jupyterlab

Developer installation
-------------------------

First fork the git repository (https://github.com/openalea/hydroroot) and clone it locally see https://docs.github.com/en/get-started/quickstart/fork-a-repo.

Create a conda environment with the proper dependencies and activate the environment

::

    conda create -n hydroroot -c conda-forge -c openalea3 openalea.deploy openalea.mtg openalea.plantgl pandas matplotlib numpy scipy yaml pyyaml rsml
    conda activate hydroroot

That creates a conda environment called *hydroroot*, install all the dependencies and activates the environment.

They open an Ipython session, source to the src directory of your cloned project and enjoy.