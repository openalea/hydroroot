import multiprocessing

import numpy as np
import pylab
from openalea.core import alea
from datetime import datetime
import cPickle as pickle

pm = alea.load_package_manager()
df =  pm['hydroroot.demo']['compute flux'].instantiate()

def surface(n=1500, branching_delay=20, branching_variability=0.25, nude_tip_length=200):
    df.node(62).set_input(0,n)
    df.node(47).set_input(0,branching_delay)

    df.eval_as_expression(vtx_id=31)
    surface = df.node(31).output(0)
    return surface


cpu = multiprocessing.cpu_count()
pool = multiprocessing.Pool(cpu)

x = range(200,1500,100)

y = pool.map(surface,x)
#y = [surface(n) for n in x]

pylab.plot(x,y)
pylab.show()

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

params = {}
params['branching_delay'] = 20
params['branching_variability']=0.25
params['nude_tip_length']=200
simul = Simulation('surface', params,x,y)
fn = simul.dump() # compute the date, a filename and pickle in a file
simul = simul.load(fn)

x, y = simul.x, simul.y
pylab.plot(x,y)
pylab.show()

