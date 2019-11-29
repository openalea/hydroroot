"""Simulation of hydro on Col.

- Data Col
- Simulate several RSA for length in range(8,14, 0.5)
"""

import datetime
import pandas
import pickle
import numpy as np
from hydroroot.main import hydroroot_mtg, hydroroot_flow
from hydroroot import radius, analysis

from data import *
import random

##############################################################################
# Data
##############################################################################

TYPE='high'

# 1. Simulate several architecture.

mtgs = []
#for genotype in ['col',]:# 'esk11', 'esk15']:#, 'irx34', 'pip2122']:

def dump_mtg(g, length, type=TYPE):
    print 'length is : ', length
    pass


results = {}
def init():
    global results
    results['index']=[]
    results['length']=[]
    results['axfold']=[]
    results['radfold']=[]
    results['surface']=[]
    results['volume']=[]
    results['Jv']=[]
    results['i0']=[]
    results['i1']=[]


def add(index, length, axfold, radfold, surface, volume, Jv, i0, i1):
    global results
    results['index'].append(index)
    results['length'].append(length)
    results['axfold'].append(axfold)
    results['radfold'].append(radfold)
    results['surface'].append(surface)
    results['volume'].append(volume)
    results['Jv'].append(Jv)
    results['i0'].append(i0)
    results['i1'].append(i1)


def save(name='bench_%s'%TYPE):
    now = datetime.datetime.now()
    date = now.strftime("%Y_%m_%d_%H-%M")
    name = name+'_%s.txt'%date
    df = pandas.DataFrame(results, columns=['index', 'length', 'axfold', 'radfold', 'surface', 'volume', 'Jv', 'i0', 'i1'])
    df.to_csv(name, index=False)


length_values = np.arange(8, 15.5, 0.1).tolist()
length_values = (12.,)
for nb_time in range(30):
    for length in length_values:
        print 'Simulation ', length
        length = length/100.
        length_data = primary_length_law(length, type=TYPE)

        seed = random.randint(1,10000)
        g, surface, volume = hydroroot_mtg(
        seed=seed,
        primary_length=length,
        length_data=length_data,
        nude_length=0.02,
        )
        dump_mtg(g, length, type=TYPE)
        mtgs.append((g, length))

# 2. Save them for analysis
# 3. Compute k0 such that, Jv(k0/10) = 1/3 Jv(k0)
count = 0
init()

for g, length in mtgs:
    for k0 in (400,):#range(50, 1500, 100):
        g, Keq, Jv_global = hydroroot_flow(g, k0=k0,
                                           axial_conductivity_data=axial_law(genotype='col'),
                                           radial_conductivity_data=radial(k0))
        g, surface = radius.compute_surface(g)
        g, volume = radius.compute_volume(g)

        i0, i1 = analysis.intercept(g, (4.5/100, 9.5/100))

        add(count, length, 1., k0, surface, volume, Jv_global, i0, i1)
        print 'Simu, ', count
        count += 1
save()


