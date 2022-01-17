###############################################################################
# Authors: C. Pradal, Y. Boursiac
# Date : 14/10/2016
###############################################################################

from random import _hexlify, _urandom
from openalea.deploy.shared_data import shared_data
import pandas
import pylab
import numpy as np
from openalea.plantgl.all import Viewer

import hydroroot
from hydroroot import display
from hydroroot.law import histo_relative_law
from hydroroot import radius, markov, flux, conductance, measured_root
from hydroroot.main import hydroroot_flow
from hydroroot.analysis import intercept

#%gui qt
#%matplotlib qt

###############################################################################
# data
###############################################################################

share = shared_data(hydroroot, share_path='share/plant0916')
data = share.glob('length*.csv')

names = _names = ('LR_length_mm', 'relative_distance_to_tip')

df = None

length_data = []
for d in data:
    pd = pandas.read_csv(d, sep=';', header=1,
                     names=names)
    pd.sort('relative_distance_to_tip', inplace=True)
    #pd.plot(x='relative_distance_to_tip', y=['LR_length_mm'], xticks=range(0,100,10), yticks=range(0,160,10))
    length_data.append(pd)


def length_law(pd, scale_x=1/100., scale_y=1., scale=1e-4):
    """
    scale
    """
    x = pd.relative_distance_to_tip.tolist()
    y = pd.LR_length_mm.tolist()

    size = 5.*scale_x
    #pylab.clf()
    _length_law = histo_relative_law(x, y, size=size, scale_x = scale_x, scale_y=1.e-3*scale_y, scale=scale, plot=False)
    return _length_law

###############################################################################
# parameters
###############################################################################

TYPE='default_13'

primary_length = 0.13 # 15 cm
delta = 0.002
beta = 0.25 # 25 %
order_max = 4
segment_length = 1e-4
nude_length = 0.02
seed = 2

ref_radius = 1e-4 # in m
order_decrease_factor = 0.7

# parameters
k0 = 400.
Jv = 0.1
psi_e = 0.4
psi_base = 0.

# laws
acol = axial_conductivity_data = (
    [0., 0.03,  0.06, 0.09, 0.12, 0.15, 0.18],
    [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4,292.6e-4,262.5e-4,511.1e-4]
)
def radial(v=300, scale=1):
    xr = [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16]
    yr = [v*scale]*len(xr)
    return xr, yr

def axial(scale=1):
    x, y = acol
    y = [a*scale for a in y]
    return x, y

rcol = radial_conductivity_data = radial(300)

axfold=1
radfold=1


###############################################################################
# Main simulation function
###############################################################################

def my_seed():
    """ Define my own seed function to capture the seed value. """
    return int(int(_hexlify(_urandom(2500)), 16)%100000000)

def my_run(primary_length, axfold=1., radfold=1.):
    """ Simulate Arabidopsis architecture

    """
    n=None

    _seed = my_seed()
    # Report the code into a function
    nb_nude_vertices = int(nude_length / segment_length)
    branching_delay = int(delta / segment_length)

    if n is None:
        nb_vertices = int(primary_length / segment_length)
    else:
        nb_vertices = n
        warn("Use primary_length instead")

    length_max_secondary = length_data[0].LR_length_mm.max()*1e-3 # in m
    law_order1 = length_law(length_data[0], scale_x=primary_length/100., scale=segment_length)
    law_order2 = length_law(length_data[1], scale_x=length_max_secondary/100., scale=segment_length)

    g = markov.markov_binary_tree(
            nb_vertices=nb_vertices,
            branching_variability=beta,
            branching_delay=branching_delay,
            length_law=[law_order1, law_order2],
            nude_tip_length=nb_nude_vertices,
            order_max=order_max,
            seed=_seed)


    # compute radius property on MTG
    g = radius.ordered_radius(g, ref_radius=ref_radius, order_decrease_factor=order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    _length = g.nb_vertices(scale=1)*segment_length
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    i1, = intercept(g, [0.045])

    # compute axial & radial
    g, Keq, Jv_global = hydroroot_flow(g,
                                       segment_length=1.e-4,
                                       k0=k0,
                                       Jv=Jv,
                                       psi_e=psi_e,
                                       psi_base=psi_base,
                                       axial_conductivity_data=axial(axfold),
                                       radial_conductivity_data=radial(300., radfold))
    return g, axfold, radfold, _length, surface, Jv_global, i1, _seed


# for i in range(3):
#     g, Jv_global, _length, surface, axfold, radfold, i1, seed= my_run()
#     scene = display.plot(g)
#     Viewer.display(scene)
#     print Jv_global, _length, surface, i1



###############################################################################
# RUN IT!!!!!
###############################################################################

import datetime

results = {}
def init():
    global results
    results['index']=[]
    results['primary_length']=[]
    results['axfold']=[]
    results['radfold']=[]
    results['length']=[]
    results['surface']=[]
    results['Jv']=[]
    results['i0']=[]
    results['seed']=[]


def add(index, primary_length, axfold, radfold, length, surface, Jv, i0, seed):
    global results
    results['index'].append(index)
    results['primary_length'].append(primary_length)
    results['axfold'].append(axfold)
    results['radfold'].append(radfold)
    results['length'].append(length)
    results['surface'].append(surface)
    results['Jv'].append(Jv)
    results['i0'].append(i0)
    results['seed'].append(seed)


def save(name='bench_%s'%TYPE):
    now = datetime.datetime.now()
    date = now.strftime("%Y_%m_%d_%H-%M")
    name = 'results/'+name+'_%s.txt'%date
    df = pandas.DataFrame(results, columns=['index', 'primary_length', 'axfold', 'radfold', 'length', 'surface', 'Jv', 'i0', 'seed'])
    df.to_csv(name, index=False)

count = 0
init()

length_values = np.arange(12, 15.5, 0.5).tolist()
length_values = (13.,)
for nb_time in range(10000):
    for length in length_values:
        count += 1

        length = length/100.
        g, axfold, radfold, _length, surface, Jv, i1, seed= my_run(primary_length=length)

        add(count, primary_length, axfold, radfold, _length, surface, Jv, i1, seed)
        print('Simu, ', count)

# 2. Save them for analysis
# 3. Compute k0 such that, Jv(k0/10) = 1/3 Jv(k0)


# for g, length in mtgs:
#     for k0 in (400,):#range(50, 1500, 100):
#         g, Keq, Jv_global = hydroroot_flow(g, k0=k0,
#                                            axial_conductivity_data=axial_law(genotype='col'),
#                                            radial_conductivity_data=radial(k0))
#         g, surface = radius.compute_surface(g)
#         g, volume = radius.compute_volume(g)

#         i0, i1 = analysis.intercept(g, (4.5/100, 9.5/100))

#         add(count, length, 1., k0, surface, volume, Jv_global, i0, i1)
#         print 'Simu, ', count
#         count += 1
save()


