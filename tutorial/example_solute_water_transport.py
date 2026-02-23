#!/usr/bin/env python
# coding: utf-8

# # Solute and water transport parameters

# In[2]:


import math
from openalea.hydroroot import flux
from openalea.hydroroot.main import root_builder
from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.read_file import read_archi_data
from openalea.hydroroot.conductance import set_conductances
from openalea.hydroroot.water_solute_transport import pressure_calculation_no_non_permeating_solutes, init_some_MTG_properties
from openalea.widgets.plantgl import PlantGL # notebook viewer 3D
from openalea.plantgl.algo.view import view # 2D view
from openalea.hydroroot.display import mtg_scene


# Read the yaml file and set the Parameters variables, assuming that the code is run from the example folder

# In[3]:


parameter = Parameters()
parameter.read_file('parameters_Ctr-3P2.yml')


# In the code the concentration are in $mol.\mu L^{-1}$, then we convert them from $mol.m^{-3}$

# In[4]:


Cse = parameter.solute['Cse'] * 1e-9 # mol/m3 -> mol/microL, external permeating solute concentration
Ce = parameter.solute['Ce'] * 1e-9 # mol/m3 -> mol/microL, external non-permeating solute concentration


# Read the architecture file and build the MTG

# In[5]:


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


# Set the conductance in the MTG (in previous examples that was done in hydroroot_flow), set some other properties in *init_some_MTG_properties* and perform some initialization. Note that here *parameter.hydro['k0']* is a float.

# In[6]:


g = set_conductances(g, axial_pr = parameter.hydro['axial_conductance_data'], k0_pr = parameter.hydro['k0']) 
g = flux.flux(g, psi_e = parameter.exp['psi_e'], psi_base = parameter.exp['psi_base'])  # initialization
g = init_some_MTG_properties(g, tau = parameter.solute['J_s'], Cini = parameter.solute['Cse'])


# Perform the calculation, this is a Newtown-Raphson loop on a matrix system. *pressure_calculation_no_non_permeating_solutes*, as its name indicates, is a solving function where no non-permeating solute is considered inside the root.

# In[7]:


eps = 1.0e-9 # global: stop criterion for the Newton-Raphson loop in Jv_P_calculation and Jv_cnf_calculation
nb_v = g.nb_vertices()
Fdx = 1.0
Fdx_old = 1.
while Fdx > eps:
    g, dx, data, row, col = pressure_calculation_no_non_permeating_solutes(g, sigma = parameter.solute['Sigma'], 
                                                                           Ce = Ce,
                                                                           Cse = Cse, 
                                                                           Pe = parameter.exp['psi_e'], 
                                                                           Pbase = parameter.exp['psi_base'])
    Fdx = math.sqrt(sum(dx ** 2.0)) / nb_v
    if abs(Fdx - Fdx_old) < eps: break
    Fdx_old = Fdx
Jv = g.property('J_out')[1]


# In[8]:


result=f"""
primary length (m): {primary_length}
surface (m2): {surface}
total length (m): {total_length}
flux (microL/s): {Jv}
"""
print(result)


# Display the concentration in the architecture 3D view

# In[9]:


s = mtg_scene(g, prop_cmap = 'C') # create a scene from the mtg with the property j is the radial flux in ul/s
PlantGL(s) 


# In[ ]:





# In[ ]:




