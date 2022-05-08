# -*- coding: utf-8 -*-
__revision__ = "$Id: $"

import sys
import os

from setuptools import setup, find_packages

# Define metainfo

name = 'OpenAlea.HydroRoot'
package = 'HydroRoot'
description= 'HydroRoot package for OpenAlea.'
long_description= 'OpenAlea.HydroRoot is a hydraulic root architecture modelling and a root architecture system generator package.'
authors= 'Christophe Pradal, Yann Boursiac, Mikael Lucas, Fabrice Bauget'
authors_email = 'christophe pradal at cirad fr, yann boursiac at inrae fr'
url = 'https://github.com/openalea/hydroroot'
license = 'Cecill-C'

packages=find_packages('src')
package_dir={'': 'src'}

# find version number in src/alinea/astk/version.py
versiond = {}
with open("src/hydroroot/version.py") as fp:
    exec(fp.read(), versiond)
version = versiond["__version__"]

# List of top level wralea packages (directories with __wralea__.py) 
#wralea_entry_points = ['%s = %s'%(pkg,namespace + '.' + pkg) for pkg in top_pkgs]

setup_requires = ['openalea.deploy']

setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
    author=authors,
    author_email=authors_email,
    url=url,
    license=license,
    keywords = 'openalea, hydraulic, plant',	

    # package installation
    packages= packages,	
    package_dir= package_dir,

    # Namespace packages creation by deploy
    #namespace_packages = [namespace],
    #create_namespaces = False,
    zip_safe= False,

    # Dependencies
    setup_requires = setup_requires,

    # Eventually include data in your package
    # (flowing is to include all versioned files other than .py)
    include_package_data = True,
    # (you can provide an exclusion dictionary named exclude_package_data to remove parasites).
    # alternatively to global inclusion, list the file to include   
    #package_data = {'' : ['*.pyd', '*.so'],},

    # postinstall_scripts = ['',],

    # Declare scripts and wralea as entry_points (extensions) of your package 
    entry_points = { 
        'wralea' : ['hydroroot = hydroroot_wralea',
                    'hydroroot.demo = hydroroot_wralea.demo',
                    'hydroroot.data= hydroroot_wralea.data',],
        },
    )


