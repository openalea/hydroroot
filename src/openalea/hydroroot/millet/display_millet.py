""" TODO: Header

"""

########################################################################################
# Imports

import random
from collections import defaultdict

import numpy as np

from openalea.mtg.mtg import *
from openalea.mtg.io import *
from openalea.mtg.traversal import *
from openalea.mtg.turtle import *
from openalea.mtg.plantframe import*
from openalea.mtg.plantframe.color import *
from openalea.plantgl import *
import openalea.plantgl.all as pgl

from openalea.hydroroot.millet import turtle, law
from openalea.hydroroot.display import my_colormap
########################################################################################


def plot_mtg(g, static, 
             save_image=False, 
             dir='C:/Users/ndour/Desktop/mon_dossier/my_Ph.D/Plant_images_BD/Film_MTG/architecture'):
    """
    Plot an MTG at final stage or during its development.

    Parameters:
    ===========
      - g : MTG
      - state : static or dynamic

    Example:
    ========

      >>> plot_mtg(g, "static")
      >>> plot_mtg(g, "dynamic")

    """
    g =law.developmental_age(g,15)
    ages = g.property('age')
    age_max = max(ages.values())
    l= np.linspace(1,age_max,50)
    l=map(int,l)
    if static == True:
        g, scene = turtle.mtg_turtle_time(g,age_max)
        Viewer.display(scene)
    else:
        for t in l:
            g, scene = turtle.mtg_turtle_time(g,t)
            Viewer.display(scene)
            if save_image==True:
                Viewer.saveSnapshot(dir+'/millet_%04d.png'%t)

#############################################################################################################

def plot_property(g, static=True,
                  save_image=False,
                  prop_cmap='radius', 
                  cmap='jet',
                  lognorm=True,
                  dir='C:/Users/ndour/Desktop/mon_dossier/my_Ph.D/Plant_images_BD/Film_MTG/architecture+j'):

    #r_base, r_tip = float(r_base), float(r_tip)

    #if not has_radius:
        #radius.discont_radius(g,r_base=r_base, r_tip=r_tip)

    #turtle = turt.PglTurtle()
    #turtle.down(180)
    #scene = turt.TurtleFrame(g, visitor= visitor, turtle=turtle, gc=False)

    age = g.property('age')
    age_max = max(age.values())

    l= np.linspace(1,age_max,75)
    l=map(int,l)

    if static == True:
        g, scene = turtle.mtg_turtle_time(g,age_max)

        # Compute color from radius
        my_colormap(g,prop_cmap, cmap=cmap, lognorm=lognorm)

        # F. Bauget 2025-08-28: WIP python 2 to 3, got "AttributeError: 'Shape' object has no attribute 'getId'"
        # shapes = dict( (sh.getId(),sh) for sh in scene)
        shapes = scene.todict()

        colors = g.property('color')
        for vid in shapes:
            # shapes[vid].appearance = pgl.Material(colors[vid])
            for sh in shapes[vid]:
                sh.appearance = pgl.Material('Color%d' % vid, colors[vid])
        # scene = pgl.Scene(shapes.values())
        scene = pgl.Scene([sh for shid in shapes.values() for sh in shid])
        Viewer.display(scene)

    else:
        for t in l:
            g, scene = turtle.mtg_turtle_time(g,t)
            # Compute color from radius
            my_colormap(g,prop_cmap, cmap=cmap, lognorm=lognorm)
            # F. Bauget 2025-08-28: WIP python 2 to 3, got "AttributeError: 'Shape' object has no attribute 'getId'"
            # shapes = dict( (sh.getId(),sh) for sh in scene)
            shapes = scene.todict()
            colors = g.property('color')
            for vid in shapes:
                # shapes[vid].appearance = pgl.Material(colors[vid])
                for sh in shapes[vid]:
                    sh.appearance = pgl.Material('Color%d'%vid, colors[vid])
            # scene = pgl.Scene(shapes.values())
            scene = pgl.Scene([sh for shid in shapes.values() for sh in shid])
            Viewer.display(scene)
            if save_image==True:
                Viewer.saveSnapshot(dir+'/millet_%04d.png'%t)



###############################################################################



