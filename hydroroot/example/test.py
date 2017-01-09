# measured laws
import numpy as np

length_data = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
axial_conductivity_data = (
    [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15, 0.16],
    [1e-4, 1e-4, 1e-4, 4e-4,4e-4,2e-3,2e-3,2.5e-3,2.5e-3,2.5e-3,10e-3, 1e-1]
)
axial_conductivity_data = (
    [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12, 0.135, 0.15, 0.16],
    #[1e-4, 1e-4, 1e-4, 3e-4, 10e-4, 1.1e-3, 1.2e-3, 1.5e-3, 2e-3, 2.4e-3, 3e-3, 1e-1],
    [1e-4, 1e-4, 1e-4, 8e-4, 28e-4, 2.9e-3, 3.2e-3, 3.5e-3, 3.8e-3, 4e-3, 10e-3,1e-1])

radial_conductivity_data = (
    [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16],
    [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300])
xr, yr = radial_conductivity_data

def f(v=300):
    y = list(np.array(yr)*(v/300.))
    return xr, y

length_data = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
primary_length = 0.1

from hydroroot.main import hydroroot
from hydroroot import display

seed = 2
g, surface, volume, Keq, Jv_global = hydroroot(
    seed=seed,
    primary_length=primary_length,
    length_data=length_data,
    axial_conductivity_data=axial_conductivity_data,
    radial_conductivity_data=radial_conductivity_data,
    )

