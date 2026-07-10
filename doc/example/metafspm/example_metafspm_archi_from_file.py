#!/usr/bin/env python
# coding: utf-8

# # Measured architecture

# Read architecture from text file and run a simulation of the sap flux from an Arabidopsis de-topped root plunged in a hydroponic solution at a hydrostatic pressure of 0.4 Mpa when its  base is at the atmospheric pressure.

# ## The text format architecture file
# 
# |distance_from_base_(mm)  |lateral_root_length_(mm)  |order|
# |:-:  | :-:  | :-: |
# |0.89                     |90.81             	       |1|
# |3.02                     |63.98             	       |1|
# |102.94                   |0.0             	         |1|
# |2.14                     |23.72             	       |1-1|
# |90.81                    |0.0             	         |1-1|
# |2.48                     |5.15             	       |1-2|
# |63.98                    |0.0             	         |1-2|
# 
# This is a tab separated text file with 3 columns:
# 
#     1. the distance from base of the branching laterals in mm
# 
#     2. the lateral root length in mm
# 
#     3. a string of one or more number indicating the parent root
# 
# In the example above, the root has two lateral of 1st order and on each of them one lateral of 2d order. The order of 1 indicates that the laterals are on the primary root. The last line with order 1 with 0.0 in the second column indicates the primary root tip. The line with order 1-1 indicates that this a second order lateral on the first lateral positioned at 2.14 mm from the branching on the primary root. And so on.
# 
# A 4th column with the averaged diameter of the root may be given, that may be used to build the MTG representing the architecture.

# ## Running the calculation
# If the package HydroRoot is not installed, the following examples can be run by cloning the sources from git and then sourcing the src directory in Ipython console for instance like this:

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


import matplotlib
from openalea.widgets.plantgl import PlantGL # notebook viewer 3D
from openalea.plantgl.algo.view import view # 2D view
from openalea.hydroroot.display import mtg_scene
from openalea.hydroroot.read_file import read_archi_data
from openalea.hydroroot.main import hydroroot_flow, root_builder
from openalea.hydroroot.generator import measured_root
from openalea.hydroroot import model
from openalea.mtg import MTG


# Read the architecture file as DataFrame

# In[3]:


df = read_archi_data('../data/plant-01.txt')


# Building the MTG from the file, and  return the primary root length, the total length and the surface. The seed refer to the seed of the root generator when the MTG is not built from a file but is generated.

# In[4]:

segment_length = 1.0e-4
time_step = 1.0e-4

g = measured_root.mtg_from_aqua_data(df, segment_length)
hydromodel = model.HydroRootModel(g, time_step)
hydromodel.compute_metrics()

# Some conductance data versus distance to tip

# In[5]:


hydromodel.compute_radial_conductance(data = ([0, 0.2],[30.0,30.0]))
hydromodel.compute_axial_conductance(data = ([0, 0.2],[3.0e-7,4.0e-4]))


# Flux and equivalent conductance calculation, for a root in an external hydroponic medium at 0.4 MPa, its base at 0.1 MPa, and with the conductances set above.

# In[6]:

# g = hydromodel.g
# g, keq, jv = hydroroot_flow(g, psi_e = 0.4, psi_base = 0.1, axial_conductivity_data = K_axial_data, radial_conductivity_data = k_radial_data)
hydromodel.psi_e = 0.4
hydromodel.psi_base = 0.1
hydromodel.hydrostatic_solver_flux()

# In[7]:


print('equivalent root conductance (microL/s/MPa): ',hydromodel.keq, 'sap flux (microL/s): ', hydromodel.Jv)


# ## Plots

# Display the local water uptake heatmap in 3D

# In[8]:


s = mtg_scene(g, prop_cmap = 'j') # create a scene from the mtg with the property j is the radial flux in ul/s
PlantGL(s) # display the root into the plantgl oawidget


# You may change the property to display to the hydrostatic pressure inside the xylem vessels for instance
# 
# to reduce notebook size we use here a 2D view but you can use the openalea.widgets `PlantGL(s)` to display an interactive 3D view

# In[9]:


s = mtg_scene(g, prop_cmap='psi_in')
view(s) #PlantGL(s)


# You may change the radial conductivity and see the impact on the water uptake

# In[10]:


hydromodel.compute_radial_conductance(data = ([0, 0.2],[300.0,300.0]))
hydromodel.hydrostatic_solver_flux()
print('sap flux (microL/s): ', hydromodel.Jv)
s = mtg_scene(g, prop_cmap='j')
view(s) # to reduce notebook size we use here a 2D view but use PlantGL(s) to 3D


# Or the axial conductance

# In[11]:


hydromodel.compute_radial_conductance(data = ([0, 0.2],[30.0,30.0]))
hydromodel.compute_axial_conductance(data = ([0, 0.2],[3.0e-7,4.0e-4]))
hydromodel.hydrostatic_solver_flux()
print('sap flux (microL/s): ', hydromodel.Jv)
s = mtg_scene(g, prop_cmap='j')
view(s) # to reduce notebook size we use here a 2D view but use PlantGL(s) to 3D


# ## PlantGl viewer
# 
# This is only working on **local machine**.
# 
# HydroRoot has also a function using the standalone Plantgl viewer allowing to specify the mtg property.
# 
# The following will open a standalone window with a 3D representation of the root. Uncomment the lines to run the cell locally.

# In[12]:


# from openalea.hydroroot.display import plot
# %gui qt
# plot(g, prop_cmap='psi_in')
from openalea.hydroroot.init_parameter import Parameters
parameter = Parameters()
parameter.read_file('parameters_plant_01.yml')
scenario = parameter.metafspm_scenario()
hydromodel = model.HydroRootModel(g, time_step, **scenario)
hydromodel.hydrostatic_solver_flux()
print('sap flux (microL/s): ', hydromodel.Jv)
