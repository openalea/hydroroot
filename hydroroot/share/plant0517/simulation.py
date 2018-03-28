###############################################################################
#
# Authors: C. Pradal, Y. Boursiac
# Date : 14/10/2016
#
###############################################################################

######
# Imports

# VERSION = 2

from random import _hexlify, _urandom

import numpy as np
import pandas
import pylab

from openalea.deploy.shared_data import shared_data
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

# column names
names = _names = ('LR_length_mm', 'relative_distance_to_tip')

# data frame
df = None

# contains 2 dataframe : order1 and order2
length_data = []
# order1 & order2
for d in data:
    pd = pandas.read_csv(d, sep=';', header=1,
                     names=names)
    pd.sort_values(by='relative_distance_to_tip', inplace=True)
    #pd.plot(x='relative_distance_to_tip', y=['LR_length_mm'], xticks=range(0,100,10), yticks=range(0,160,10))
    length_data.append(pd)


def length_law(pd, scale_x=1/100., scale_y=1., scale=1e-4, uniform=True):
    """
    scale
    """
    x = pd.relative_distance_to_tip.tolist()
    y = pd.LR_length_mm.tolist()

    # size of the windows: 5%
    size = 5.*scale_x
    #pylab.clf()
    _length_law = histo_relative_law(x, y, size=size, scale_x=scale_x, scale_y=1.e-3*scale_y, scale=scale, plot=False, uniform=uniform)
    return _length_law

###############################################################################
# parameters
###############################################################################

TYPE='default'

primary_length = 0.13 # 15 cm
delta = 0.002
beta = 0.25 # 25 %
order_max = 4
segment_length = 1e-4
nude_length = 0.02
#seed = 2

ref_radius = 7e-5 # in m
order_decrease_factor = 0.7

# parameters
#k0 = 400.
k0 = 92.

Jv = 0.1
psi_e = 0.4
psi_base = .101325

# laws
acol = axial_conductivity_data = (
    [0., 0.03,  0.06, 0.09, 0.12, 0.15, 0.18],
    [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4,292.6e-4,262.5e-4,511.1e-4]
)
def radial(v=92, scale=1):
    xr = [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16]
    yr = [v*scale]*len(xr)
    return xr, yr

def axial(scale=1):
    x, y = acol
    y = [a*scale for a in y]
    return x, y

rcol = radial_conductivity_data = radial(92)

axfold=1
radfold=1

def plot(g, **kwds):

    from IPython import get_ipython

    # ipython = get_ipython()
    ipython = None
    if ipython:
        ipython.magic('gui qt')

    from openalea.plantgl.all import Viewer
    from hydroroot.display import plot as _plot

    Viewer.display(_plot(g, **kwds))

###############################################################################
# Main simulation function
###############################################################################

ONE_LAW=False
EXPOVARIATE = True

def my_seed():
    """ Define my own seed function to capture the seed value. """
    return int(long(_hexlify(_urandom(2500)), 16)%100000000)

def my_run(primary_length, k0=k0,
           axfold=1., radfold=1.,
           seed=None,
           ref_radius=ref_radius,
           order_decrease_factor=0.7,
           delta=delta,
           nude_length=nude_length):
    """ Simulate Arabidopsis architecture

    Add parameters like delta
    """
    n=None

    if seed is None:
        _seed = my_seed()
    else:
        _seed = seed

    # Report the code into a function
    nb_nude_vertices = int(nude_length / segment_length)
    branching_delay = int(delta / segment_length)

    if n is None:
        nb_vertices = int(primary_length / segment_length)
    else:
        nb_vertices = n
        warn("Use primary_length instead")

    length_max_secondary = length_data[0].LR_length_mm.max()*1e-3 # in m

    # Just use the same order1 law
    law_order1 = length_law(length_data[0], scale_x=primary_length/100., scale=segment_length)
    if EXPOVARIATE:
        law_order2 = length_law(length_data[1], scale_x=length_max_secondary/100., scale=segment_length, uniform='expo')
    else:
        law_order2 = length_law(length_data[1], scale_x=length_max_secondary/100., scale=segment_length)

    g = markov.markov_binary_tree(
            nb_vertices=nb_vertices,
            branching_variability=beta,
            branching_delay=branching_delay,
            length_law=[law_order1, law_order2] if not ONE_LAW else law_order1,
            nude_tip_length=nb_nude_vertices,
            order_max=order_max,
            seed=_seed)


    # compute radius property on MTG
    g = radius.ordered_radius(g, ref_radius=ref_radius, order_decrease_factor=order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # _length is the total length of the RSA (sum of the length of all the segments)
    _length = g.nb_vertices(scale=1)*segment_length
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    # Compute the intercept at 1,2,3,4.5,6 and 8 cm
    # TODO : np.arange(3, 11, 1)*1e-2
    i1,i2,i3,i4,i6,i8 = intercept(g, [0.01, 0.02, 0.03, 0.045, 0.06, 0.08])

    # compute axial & radial
    g, Keq, Jv_global = hydroroot_flow(g,
                                       segment_length=1.e-4,
                                       k0=k0,
                                       Jv=Jv,
                                       psi_e=psi_e,
                                       psi_base=psi_base,
                                       axial_conductivity_data=axial(axfold),
                                       radial_conductivity_data=radial(k0, radfold))
    return g, axfold, radfold, _length, surface, Jv_global, i1,i2,i3,i4,i6,i8, _seed


# for i in range(3):
#     g, Jv_global, _length, surface, axfold, radfold, i1, seed= my_run()
#     scene = display.plot(g)
#     Viewer.display(scene)
#     print Jv_global, _length, surface, i1



###############################################################################
# RUN IT!!!!!
# Startegies for saving intermediate output of the model
###############################################################################

import datetime

results = {}
def init():
    global results
    results['index'] = []
    results['primary_length'] = []
    results['k0'] = []
    results['axfold'] = []
    results['radfold'] = []
    results['length'] = []
    results['surface'] = []
    results['Jv'] = []
    results['i1'] = []
    results['i2'] = []
    results['i3'] = []
    results['i4'] = []
    results['i6'] = []
    results['i8'] = []
    results['delta'] = []
    results['nude_length'] = []
    results['seed'] = []


def add(index, primary_length, k0, axfold, radfold, delta, nude_length, length, surface, Jv, i1,i2,i3, i4, i6, i8, seed):
    global results
    results['index'].append(index)
    results['primary_length'].append(primary_length)
    results['k0'].append(k0)
    results['axfold'].append(axfold)
    results['radfold'].append(radfold)
    results['delta'].append(delta)
    results['nude_length'].append(nude_length)
    results['length'].append(length)
    results['surface'].append(surface)
    results['Jv'].append(Jv)
    results['i1'].append(i1)
    results['i2'].append(i2)
    results['i3'].append(i3)
    results['i4'].append(i4)
    results['i6'].append(i6)
    results['i8'].append(i8)
    results['seed'].append(seed)


def save(name='bench_%s'%TYPE, xls=False):
    now = datetime.datetime.now()
    date = now.strftime("%Y_%m_%d_%H-%M")
    df = pandas.DataFrame(results, columns=['index', 'primary_length', 'k0', 'axfold', 'radfold', 'delta', 'nude_length', 'length', 'surface',
                                            'Jv', 'i1','i2','i3', 'i4','i6','i8', 'seed'])
    if not xls:
        name = name+'_%s.csv'%date
        df.to_csv(name, index=False)
    else:
        name = name+'_%s.xlsx'%date
        writer = pandas.ExcelWriter(name)
        df.to_excel(writer,'Sheet1')

def main():
    count = 0
    init()

    # values are from 10 to 15 cm
    length_values = np.arange(10, 15, .5).tolist()
    deltas = np.arange(1.e-3, 2.5e-3, 1e-4).tolist()
    nudes = np.arange(0.009, 0.038, 3.e-3 ).tolist()

    axfolds = [1./16, 1./8, 1./4, 1./2, 1., 2., 3., 4., 5.]
    k0 = 92
    radfolds = [1./4, 1/2., 1., 2., 4., 8.]

    axfolds = [1.]
    radfolds = [1.]
    #length_values = (13.,)
    N = 1
    # run 200 times the model

    # predict the number of simulation run
    nb_steps = N * len(length_values) * len(deltas) * len(axfolds) * len(radfolds) * len(nudes)
    print 'Simulation runs: ', nb_steps
    print '#############################'

    #to comment
    #return


    for nb_time in range(N):
        for length in length_values:
            for delta in deltas:
                for axfold in axfolds:
                    for radfold in radfolds:
                        for nude_length in nudes:
                            count += 1

                            # convert to meters
                            plength = length/100.
                            g, axfold, radfold, _length, surface, Jv, i1, i2, i3,i4,i6,i8, seed = my_run(primary_length=plength,
                                k0=k0, axfold=axfold, radfold=radfold, delta=delta, nude_length=nude_length)

                            add(count, plength, k0, axfold, radfold, delta, nude_length, _length, surface, Jv,
                                i1,i2,i3,i4,i6,i8, seed)
                            print 'Simu, ', count

    save(xls=False)


def reproduce(fn='171002_input_seeds_individuals.txt'):
    rep_names = ('primary_length', 'k0', 'axfold', 'radfold', 'delta', 'nude_length', 'seed')
    filename = fn
    rep_pd = pandas.read_csv(filename, sep=';', header=0,
                             names=rep_names)

    r_length = rep_pd.primary_length.tolist()
    r_k0 = rep_pd.k0.tolist()
    r_axfold = rep_pd.axfold.tolist()
    r_radfold = rep_pd.radfold.tolist()
    r_delta = rep_pd.delta.tolist()
    r_nude = rep_pd.nude_length.tolist()
    r_seed = rep_pd.seed.tolist()

    number_of_runs = len(r_length)

    init()

    for i in range(number_of_runs):
        length = r_length[i]
        k0 = r_k0[i]
        axfold = r_axfold[i]
        radfold = r_radfold[i]
        delta = r_delta[i]
        nude_length = r_nude[i]
        seed = int(r_seed[i])
        print "length, seed ", length, seed
        g, axfold, radfold, _length, surface, Jv, i1, i2, i3, i4, i6, i8, seed = my_run(primary_length=length,
                                                                    k0=k0,
                                                                    axfold=axfold,
                                                                    radfold=radfold,
                                                                    delta=delta,
                                                                    nude_length=nude_length,
                                                                    seed=seed)


        add(i, length, k0, axfold, radfold, delta, nude_length, _length, surface, Jv, i1, i2,
            i3,i4,i6,i8, seed)
        #plot(g)
        #raw_input("Press Enter to continue...")

    save('reproduce_bench_%s'%TYPE)



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

#main()
#reproduce()