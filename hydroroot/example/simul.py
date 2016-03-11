from openalea.deploy.shared_data import shared_data
import pandas
import pylab
import hydroroot

share = shared_data(hydroroot, share_path='share')
data = share.glob('length*.csv')[1:]

def read_data(data=data):
    names = _names = ('relative_position', 'internode_length', 'LR_length', 'distance_to_tip')

    df = None
    for d in data:
        pd = pandas.read_csv(d, sep=';', header=1,
                         names=names)
        pd.sort('distance_to_tip', inplace=True)
        pd.LR_length.cumsum()
        pd['cum'] = pd.LR_length.cumsum()
        pd['length'] = pd.internode_length.cumsum()



        #pd.plot(x='length', y=['internode_length'],xticks=range(0,135,10), yticks=range(0,4))
        if df is None:
            df = pd
        else:
            df = df.append(pd)

        pd = df


    pd['cumsum'] = pd.LR_length.cumsum()
    pd.sort('distance_to_tip', inplace=True)
    #pd.plot(x='distance_to_tip', y=['cumsum'],xticks=range(0,135,10))#, 'cum'],xticks=range(0,101,10))

    return pd


from law import multi_law

pd = read_data(data)

def length_law(intensity='mean', pd=pd):
    x = pd.distance_to_tip.tolist()
    y = pd.LR_length.tolist()
    scale_x = 0.16/100.
    size = 5.*scale_x
    #pylab.clf()
    length_min, length_max, length_mean = multi_law(x, y, size=size, scale_x = scale_x, scale_y=1.e-3, plot=True)
    if intensity == 'mean':
        return length_mean
    elif intensity == 'min':
        return length_min
    else:
        return length_max


from openalea.mtg.traversal import pre_order2


def nb_root(g, l):
    length= {}
    root = 1
    dl = 1.e-4

    if 'mylength' in g.property_names():
        length = g.property('mylength')
    else:
        for v in pre_order2(g,root):
            pid = g.parent(v)
            length[v] = length[pid]+dl if pid else dl
        g.properties()['mylength'] = length

    count = 0
    for v in g:
        pid = g.parent(v)
        if pid and (length[pid] <= l <= length[v]):
            count+=1
    return count



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

#rcol = radial_conductivity_data = radial(300)

from hydroroot.main import hydroroot


def run_simulation(primary_length, axfold, radfold):
    # Generate plant

    g, surface, volume, Keq, Jv_global = hydroroot(
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
    Jv = Jv,
    psi_e = psi_e,
    psi_base=psi_base,
    length_data=length_law('mean'),
    axial_conductivity_data=axial(axfold),
    radial_conductivity_data=radial(300., radfold))


    dist = (0.045, 0.055, 0.095, 0.105)
    intercepts = [nb_root(g, x) for x in dist]



    return primary_length, surface, volume, Jv_global, intercepts


results = {}
def init():
    global results
    results['length']=[]
    results['axfold']=[]
    results['radfold']=[]
    results['surface']=[]
    results['volume']=[]
    results['Jv']=[]
    results['i0']=[]
    results['i1']=[]
    results['i2']=[]
    results['i3']=[]


def add(length, axfold, radfold, surface, volume, Jv, i0, i1, i2, i3):
    global results
    results['length'].append(length)
    results['axfold'].append(axfold)
    results['radfold'].append(radfold)
    results['surface'].append(surface)
    results['volume'].append(volume)
    results['Jv'].append(Jv)
    results['i0'].append(i0)
    results['i1'].append(i1)
    results['i2'].append(i2)
    results['i3'].append(i3)

def save(name='bench.txt'):
    df = pandas.DataFrame(results, columns=['length', 'axfold', 'radfold', 'surface', 'volume', 'Jv', 'i0', 'i1', 'i2', 'i3'])
    df.to_csv(name, index=False)

count = 0
for i in range(5):
    for length in (0.10, 0.12, 0.13, 0.14, 0.15):
        print 'Lenght ', length
        for axfold in (1., 4., 8., 16., 32. ):
            count += 1
            init()
            for radfold in (1., 4., 8., 16., 32.):
                print 'run simu'
                length, surface, volume, jv, intercepts = run_simulation(length, axfold, radfold)
                add(length, axfold, radfold, surface, volume, jv, *intercepts)
            save('bench%d.txt'%count)

