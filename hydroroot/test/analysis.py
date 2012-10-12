import multiprocessing

import numpy as np
import pylab
from openalea.core import alea
from datetime import datetime
import cPickle as pickle

# Load a visualea graph as a script object
pm = alea.load_package_manager()
df =  pm['hydroroot.demo']['compute flux'].instantiate()

# Functions to set parameter input for visualea nodes & Get the evaluation result from a visualea node
def surface(root_length=1500, branching_delay=20, branching_variability=0.25, nude_tip_length=200):
    df.node(47).set_input(0,root_length)
    df.node(49).set_input(0,branching_delay)

    df.eval_as_expression(vtx_id=31)
    surface = df.node(31).output(0)
    return surface

def flux(root_length=1500, branching_delay=20, elementary_k=300, psi_e=0.4, radius_reduction_factor=0.5):
    df.node(58).set_input(0,True) # model inversion
    df.node(47).set_input(0,root_length)
    df.node(49).set_input(0,branching_delay)
    df.node(29).set_input(0,elementary_k)
    df.node(4).set_input(0,psi_e)
    df.node(34).set_input(0,radius_reduction_factor)

    df.eval_as_expression(vtx_id=6)  # force radial conductance computation
    df.eval_as_expression(vtx_id=21) # force axial conductance computation
    df.eval_as_expression(vtx_id=7)  # flux computation
    g = df.node(7).output(0)
    try:
        v_base = g.component_roots_at_scale_iter(g.root, scale=g.max_scale()).next()
    except:
        v_base = g.component_roots_at_scale(g.root, scale=g.max_scale()).next()

    J_out = g.property('J_out')
    return J_out[v_base]

def my_flux(d):
    return flux(**d)

###################################################################################
# Sensitivity analysis

from sensitivity import *

repeat = 5
factors = """
root_length
branching_delay
elementary_k
radius_reduction_factor
""".split()

intervals = {}
intervals['root_length'] = (200, 1500)
intervals['branching_delay'] = (15,50)
intervals['elementary_k'] = (100, 500)
intervals['radius_reduction_factor'] = (0.2, 0.8)

binf = [intervals[f][0] for f in factors ]
bsup = [intervals[f][1] for f in factors ]

m, space = Morris(repeat,factors, binf, bsup)

n = len(space['root_length'])

params = [dict(zip(factors, (space[f][i] for f in factors))) for i in range(n)]

y = []
for i,p in enumerate(params):
    print 'simulation ', i
    f=my_flux(p)
    y.append(f)


modalities, mu ,sigma = Morris_IS(m,y)
plotSens(modalities)



# Use as many processors as possible for computation
cpu = multiprocessing.cpu_count()
pool = multiprocessing.Pool(cpu)

y = pool.map(my_flux,params)









# Define the range of variation for tested parameter, compute the output for each value & show it as a curve
x = range(200,1600,100)
#y = pool.map(surface,x)
y = pool.map(flux,x)
#y = [surface(n) for n in x]
pylab.plot(x,y)
pylab.show()


# Object to store simulation parameters & results
class Simulation(object):
    def __init__(self, name, parameters, x, y):
        self.name = name
        self.parameters = parameters
        self.x = x
        self.y = y

    def dump(self):
        fn = self.name+'_'+str(datetime.date(datetime.now()))
        f = open(fn,'wb')
        pickle.dump(self,f)
        f.close()
        return fn

    @staticmethod
    def load(fn):
        f = open(fn, 'rb')
        obj = pickle.load(f)
        f.close()
        return obj

# Build the storage object, store the simulation parameters and results
params = {}
params['root_length']=range(200,1600,100)
params['branching_delay']=20
params['nude_tip_length']=200
params['elementary_k']=300
params['psi_e']=0.4
simul = Simulation('flux', params,x,y)
fn = simul.dump() # compute the date, a filename and pickle in a file

# Load a simulation and plot its results
simul = simul.load(fn)
x, y = simul.x, simul.y
pylab.plot(x,y)
pylab.show()

