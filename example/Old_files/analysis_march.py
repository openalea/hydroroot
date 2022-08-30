"""

"""
from openalea.deploy.shared_data import shared_data
import pandas
import pylab

import hydroroot
from hydroroot.main import hydroroot_from_data
from hydroroot.generator.measured_root import mtg_builder
from hydroroot.length import fit_law
from hydroroot import radius, markov, flux, conductance

from openalea.mtg.traversal import pre_order2


plant_id = 8


share = shared_data(hydroroot, share_path='share')
#data = (share/'plants').glob('length*.csv')[1:]
data = [share/'plants0216'/'length%d.csv'%plant_id]


def read_data(data=data):
    """
    absolute_position: distance from base
    internode_length
    LR_length
    distance_to_tip

    All values are in mm.
    """
    names = _names = ('absolute_position', 'internode_length', 'LR_length', 'distance_to_tip')

    df = None
    for d in data:
        pd = pandas.read_csv(d, sep=';', header=1,
                         names=names)

        for name in names:
            pd[name]/=1e3

        pd['length'] = pd.internode_length.cumsum()

        #pd.sort('distance_to_tip', inplace=True)
        pd.LR_length.cumsum()
        pd['cum'] = pd.LR_length.cumsum()




        #pd.plot(x='length', y=['internode_length'],xticks=range(0,135,10), yticks=range(0,4))
        if df is None:
            df = pd
        else:
            df = df.append(pd)

        pd = df


    pd['cumsum'] = pd.LR_length.cumsum()
    pd.sort('absolute_position', inplace=True)
    #pd.plot(x='distance_to_tip', y=['cumsum'],xticks=range(0,135,10))#, 'cum'],xticks=range(0,101,10))

    return pd

def nb_root(g, l):
    length= {}
    root = next(g.component_roots_at_scale_iter(g.root, scale=g.max_scale()))
    dl = 1.e-4

    if 'mylength' in g.property_names():
        length = g.property('mylength')
    else:
        for v in pre_order2(g,root):
            pid = g.parent(v)
            length[v] = length[pid]+dl if pid else dl
        g.properties()['mylength'] = length

    length = g.property('mylength')
    count = 0
    for v in g:
        pid = g.parent(v)
        if pid and (length[pid] <= l <= length[v]):
            count+=1
    return count


def pdata(plant_id):
    data = [share/'plants0216'/'length%d.csv'%plant_id]
    return read_data(data)



pd = read_data(data)


#primary_length = pd.length.max()
#primary_length_data = pd.length.tolist()

#lateral_length_data = pd.LR_length.tolist()

delta = 0.002
beta = 0.25 # 25 %
order_max = 4
segment_length = 1e-4
nude_length = 0.02
seed = 2

ref_radius = 1e-4 # in m
order_decrease_factor = 1.

# parameters
k0 = 300
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

length_data = ([0., 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125, 0.125],
[0., 0.001465, 0.003725, 0.01336, 0.02051, 0.023, 0.0215, 0.0399, 0.06238, 0.0945])

def run_simulation(pd, axfold, radfold, order_max=3):
    # Generate plant
    primary_length = pd.length.max()
    primary_length_data = pd.length.tolist()

    lateral_length_data = pd.LR_length.tolist()

    g, surface, volume, Keq, Jv_global = hydroroot_from_data(
    primary_length = primary_length,
    delta = delta,
    beta = beta,
    order_max = order_max,
    segment_length = segment_length,
    nude_length = nude_length,
    seed = seed,
    ref_radius = ref_radius,
    order_decrease_factor = order_decrease_factor,
    k0 = k0,
    #Jv = Jv,
    psi_e = psi_e,
    psi_base=psi_base,
    length_data=length_data,
    axial_conductivity_data=axial(axfold),
    radial_conductivity_data=radial(k0, radfold),
    primary_length_data=primary_length_data,
    lateral_length_data=lateral_length_data,
    )


    dist = (0.045, )
    intercepts = [nb_root(g, x) for x in dist]


    return g, primary_length, Jv_global, intercepts[0]


# g, length, jv, intercepts = run_simulation(pd, axfold=1, radfold=1)
#print 'length, jv, intercepts', length, jv, intercepts


results = {}
def init():
    global results
    results['plant']=[]
    results['length']=[]
    results['axfold']=[]
    results['radfold']=[]
    results['Jv']=[]
    results['intercept']=[]


def add(plant, length, axfold, radfold, Jv, i0):
    global results
    results['plant'].append(plant)
    results['length'].append(length)
    results['axfold'].append(axfold)
    results['radfold'].append(radfold)
    results['Jv'].append(Jv)
    results['intercept'].append(i0)

def save(name='bench.txt'):
    df = pandas.DataFrame(results, columns=['plant', 'length', 'axfold', 'radfold', 'Jv', 'intercept'])
    df.to_csv(name, index=False)


count = 0
for i in (1, 2, 4, 5, 7, 8):
    plant_data = pdata(i)
        #init()
    for axfold in (1., ):
        for radfold in (1., ):
            print('run simu')
            g, length, jv, intercept = run_simulation(plant_data, axfold, radfold)
            print(length, jv, intercept)
                #add(plant_id, length, axfold, radfold, jv, intercept)
#    save('bench%d.txt'%plant_id)


