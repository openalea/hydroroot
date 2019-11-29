# !!! script has to be run from hydroroot\test\ !!!

pylab qt

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
def surface(root_length=1500, branching_delay=20, branching_variability=0.25, nude_tip_length=200,radius_reduction_factor=0.7):
    df.node(47).set_input(0,root_length)
    df.node(49).set_input(0,branching_delay)
    df.node(34).set_input(0,radius_reduction_factor)

    df.eval_as_expression(vtx_id=31)
    surface = df.node(31).output(0)
    return surface

def flux(root_length=1500, branching_delay=20, elementary_k=300, psi_e=0.4, radius_reduction_factor=0.7, shuffling=False, CONSTANT=1.0, direct_input=None):
    df.node(58).set_input(0,True) # model inversion
    df.node(47).set_input(0,root_length)
    df.node(49).set_input(0,branching_delay)
    df.node(29).set_input(0,elementary_k)
    df.node(4).set_input(0,psi_e)
    df.node(34).set_input(0,radius_reduction_factor)
    df.node(63).set_input(0,shuffling)
    df.node(7).set_input(7,CONSTANT)
    df.node(21).set_input(7,direct_input)
    if direct_input is not None:
        print "axial K exp plan :"
        for key in sorted(direct_input.keys()):
            print "%s: %s" % (key, direct_input[key])

    #df.eval_as_expression(vtx_id=6)  # force radial conductance computation
    #df.eval_as_expression(vtx_id=21) # force axial conductance computation
    df.eval_as_expression(vtx_id=7)  # flux computation
    g = df.node(7).output(0)
    try:
        v_base = g.component_roots_at_scale_iter(g.root, scale=g.max_scale()).next()
    except:
        v_base = g.component_roots_at_scale(g.root, scale=g.max_scale()).next()

    J_out = g.property('J_out')
    
    if direct_input is not None :
        return J_out[v_base], direct_input
    else :
        return J_out[v_base]


def my_flux(d):
    return flux(**d)


Kref = {0.:1.00E-4,0.015:1.00E-4,0.03:1.00E-4,0.045:4.00E-4,0.06:4.00E-4,0.075:2.00E-3,0.09:2.00E-3,0.105:2.50E-3,0.135:2.50E-3,0.15:1.00E-2,0.16:1.00E-1}

def dictpanel(d, key, panel):
    dictpanel = []
    for val in panel:
        dtemp = d.copy()
        dtemp[key]=val
        dictpanel.append(dtemp)
    return dictpanel

def expplan(d):
    """
    Take a dict and return a list of dict with each key values changed by factors ranging
    from 1/1000 to 1000. Used as an generator of experiment plans for sensitivity analysis.
    """
    plan = []
    for k in d.keys():
        val = d[k]
        panel = [val/1000.,val/100.,val/10.,val,val*10.,val*100.,val*1000.]
        plan.extend(dictpanel(d,k,panel))
    return plan

plan = expplan(Kref)

def analysis():
    nb_exp = len(plan)
    output = []
    for p in plan:
        print 'iteration %s of %s' % (plan.index(p),nb_exp)
        output.append(flux(direct_input=p))
    return output

###################################################################################
# Sensitivity analysis

from sensitivity import *

repeat = 10
factors = """
root_length
branching_delay
elementary_k
radius_reduction_factor
""".split()

factors = """
elementary_k
radius_reduction_factor
""".split()

intervals = {}
#intervals['root_length'] = (200, 1500)
#intervals['branching_delay'] = (18,50)
intervals['elementary_k'] = (0.1, 500)
intervals['radius_reduction_factor'] = (0.1, 1.)

binf = [intervals[f][0] for f in factors ]
bsup = [intervals[f][1] for f in factors ]

m, space = Morris(repeat,factors, binf, bsup)

n = len(space['elementary_k'])

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

