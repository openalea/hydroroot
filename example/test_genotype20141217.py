# measured laws
import numpy as np

length_data = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
axial_conductivity_data = (
    [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15, 0.16], 
    [1e-4, 1e-4, 1e-4, 4e-4,4e-4,2e-3,2e-3,2.5e-3,2.5e-3,2.5e-3,10e-3, 1e-1]
)
axial_conductivity_data = (
    [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15, 0.16], 
    [1e-4, 1e-4, 1e-4, 3e-4, 10e-4, 1.1e-3, 1.2e-3, 1.5e-3, 2e-3, 2.4e-3, 3e-3, 1e-1], 
    [1e-4, 1e-4, 1e-4, 8e-4, 28e-4, 2.9e-3, 3.2e-3, 3.5e-3, 3.8e-3, 4e-3, 10e-3,1e-1])
    
radial_conductivity_data = (
    [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16], 
    [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300])
xr, yr = radial_conductivity_data

def f(v=300):
    y = list(np.array(yr)*(v/300.))
    return xr, y

# axial
a={}
a['col'] = (
    [0., 0.03,  0.06, 0.09, 0.12, 0.15, 0.18], 
    [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4,292.6e-4,262.5e-4,511.1e-4]
)

a['esk11'] = (
    [0., 0.03,  0.06, 0.09, 0.12], 
    [2.3e-4, 15.9e-4, 158.9e-4, 216.4e-4,192.e-4]
)

a['esk15'] = (
    [0., 0.03,  0.06, 0.09, 0.12, 0.15], 
    [0.7e-4, 3.9e-4, 168.8e-4, 220.5e-4,239.3e-4, 219.5e-4]
)

a['irx34'] = (
    [0., 0.03,  0.06, 0.09, 0.12], 
    [5.2e-4, 11.4e-4, 19.8e-4, 119.6e-4,359.9e-4]
)

a['pip2122'] = a['col']
axial_conductivity_data = a['col']


# radial
r = {}
r['col'] = f(300)
r['pip2122'] = f(239)
r['esk11'] = f(373)
r['esk15'] = f(518)
r['irx34'] = f(300)
radial_conductivity_data = r['col']

# length
l={}
nn={}
l['col'] = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
nn['col'] = 1500

l['esk11'] = [0., 0.02, 0.05, 0.09, 0.13], [0., 0., 0.01, 0.045, 0.08]
nn['esk11'] = 1200

l['esk15'] = [0., 0.02, 0.08, 0.10, 0.16], [0., 0., 0.01, 0.02, 0.05]
nn['esk15'] = 1500 

l['irx34'] = [0., 0.02, 0.05, 0.09, 0.13], [0., 0., 0.005, 0.03, 0.055]
nn['irx34'] = 1200

l['pip2122'] = l['col']
nn['pip2122'] = 1500

length_data = l['col']
n = nn['col']

from hydroroot.main import hydroroot 

for genotype in ['col', 'esk11', 'esk15', 'irx34', 'pip2122']:
    length_data = l[genotype]
    axial_conductivity_data = a[genotype]
    radial_conductivity_data = r[genotype]
    n = nn[genotype]
    seed = 2
    g, surface, volume, Keq, Jv_global = hydroroot(
    seed=seed,
    n = n,
    length_data=length_data, 
    axial_conductivity_data=axial_conductivity_data, 
    radial_conductivity_data=radial_conductivity_data,
    )
    print genotype, 'volume:', volume,'Keq:', Keq, 'Jv:', Jv_global, Jv_global/volume

"""
col : 1.02133448396 0.306400345189 1087538.56179
esk11 : 0.277413069454 0.0832239208361 2307132.08605 (-73% +112%)
esk15 : 0.163983327915 0.0491949983746 2260361.87944 (-84% +107%)

irx34 : 0.130030783485 0.0390092350454 2049689.80902 (-87% +88%)
pip2122 : 0.950280346469 0.285084103941 1011878.61326 (-7% -6%)
"""