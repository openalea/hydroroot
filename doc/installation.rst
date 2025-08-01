============
Installation
============

First install miniforge following the instructions given here https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

User installation
---------------------

HydroRoot may by installed simply on a conda environments:

::

    mamba create -n hydroroot -c conda-forge -c openalea3 openalea.hydroroot
    mamba activate hydroroot

That creates a conda environment called *hydroroot*, install in it *openalea.hydroroot* with all the dependencies and
activate the environment. Then just open an Ipython session and enjoy.

If you want notebook support, run for example:

::

    mamba install jupyterlab

Developer installation
-------------------------

First fork the git repository (https://github.com/openalea/hydroroot) and clone it locally see https://docs.github.com/en/get-started/quickstart/fork-a-repo.

Just run the following command:

::

    mamba create -f conda/environment.yml
    mamba activate hydroroot
    pip install -e .

This will first create a conda environment called *hydroroot* with the proper dependencies, then the environment will be activated,
and finally openalea.hydroroot will be installed in develop mode. As above to have notebook support run `mamba install jupyterlab`.
