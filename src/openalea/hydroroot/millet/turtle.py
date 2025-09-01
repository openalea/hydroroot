from openalea.mtg.mtg import *
from openalea.mtg.io import *
from openalea.mtg.traversal import *
from openalea.mtg.turtle import *
from openalea.mtg.plantframe import*
from openalea.mtg.plantframe.color import *
from openalea.plantgl import *


from openalea.hydroroot.millet import architecture
import random
import numpy as np

from collections import defaultdict



############################################################
# TODO : Move this code to mtg (CPL)


LR_root_angle= defaultdict(lambda: random.normalvariate(53.52674,15.14071)) #we store the lateral root angles in a default dictionary
#average_CR_angle= 20.03194444
average_CR_angle= 30.03194444
segment_length = 1e-4

tort = defaultdict(lambda: random.randint(-1,1))





def mtg_turtle_time(g, time, update_visitor=None, time_property='age' ):
    ''' Compute the geometry on each node of the MTG using Turtle geometry. 
    
    Update_visitor is a function called on each node in a pre order (parent before children).
    This function allow to update the parameters and state variables of the vertices.
    
    :Example:

        >>> def grow(node, time):

    '''

    g.properties()['geometry'] = {}
    orders = g.property('order')
    radius = g.property('radius')
    max_scale = g.max_scale()

    start_dates = g.property(time_property)
    
    #seminal = g.Trunk(g.root, Scale=1)
    
    #
    #lats= [g.Descendants(g.Successor(v)) for v in g if g.edge_type(v)=='+' and n.label !='Crown']
    
    def visitor(g, v, turtle,time):

        n = g.node(v)
        radius = n.radius * 1.e1
        length = n.length * 1.e3
        #LR_root_angle=30
        if g.edge_type(v) == '+':         #For crown root angles, the size is weak, so we use a mean to estimate their branching
            if n.label=='Crown':
                turtle.down(average_CR_angle)
                


            else: #We use here a normal law to simulate the LR angle values seeing that they follow a Gaussian distribution
                random.seed(v)
                turtle.down(LR_root_angle[v])

        #use rather curves to define crown roots tortuosity
        # if n.label=='Crown':
        #     random.seed(v)
        #     turtle.down(tort[v])


        # if n.label == 'Seminal':
        #     #Define SR tortuosity
        #     if n.order == 0:
        #         random.seed(v)
        #         turtle.down(random.random())

        #     # #Define LR tortuosity
        #     else:
        #         random.seed(v)
        #         turtle.down(random.random())

        if n.label == 'collet':
            turtle.rollL(90.)

        turtle.setId(v)
        turtle.setWidth(radius)
        #for c in n.children():
            #if c.edge_type() == '+':
                #turtle.rollL(130)
        turtle.setColor(n.order+1)
        if n.label=='Crown':
            turtle.setColor(n.order+2)
        turtle.rollL()
        turtle.F(length)
        pt = turtle.getPosition()
        n.xyz = pt

    def traverse_with_turtle_time(g, vid, time, visitor=visitor):
        turtle = PglTurtle()
        def push_turtle(v):
            try:
                start_tt = start_dates.get(v, time+1)
                if start_tt > time:
                    return False
            except: 
                pass
            if g.edge_type(v) == '+':
                turtle.push()
                turtle.setId(v)
            return True

        def pop_turtle(v):
            try:
                start_tt = start_dates.get(v, time+1)
                if start_tt > time:
                    return False
            except: 
                pass
            if g.edge_type(v) == '+':
                turtle.pop()

        if start_dates[vid] <= time:
            visitor(g,vid,turtle,time)
            #turtle.push()
        plant_id = g.complex_at_scale(vid, scale=1)
        for v in pre_order2_with_filter(g, vid, None, push_turtle, pop_turtle):
            if v == vid: continue
            # Done for the leaves
            if start_dates.get(v,time+1) > time:
                print('Do not consider ', v, time)
                continue
            visitor(g,v,turtle,time)

        scene = turtle.getScene()
        return scene, g

    for plant_id in g.component_roots_at_scale_iter(g.root, scale=max_scale):
        scene, g = traverse_with_turtle_time(g, plant_id, time)
    return g, scene


def Plot(scene):
    Viewer.display(scene)

def compute_angle(g):
    pass

#############################################################


def test(save_image= False, dir='C:/Users/ndour/Desktop/mon dossier/my_Ph.D/Plant_images_BD/Film_MTG'):
    from openalea.hydroroot.millet import architecture
    g = architecture.millet_mtg(nb_vertices=200)
    g = architecture.developmental_age(g,15)
    compute_angle(g)

    ages = g.property('age')
    age_max = max(ages.values())

    #for t in range(1, age_max, 1):
    for t in list(np.linspace(1,age_max,100)):
        g, scene = mtg_turtle_time(g,t)
        Plot(scene)
        if save_image==True:
            Viewer.saveSnapshot(dir+'/millet_%04d.png'%t)
    # sous ipython
    # %gui qt


#-------Display of the MTG
#test(save_image=False)

