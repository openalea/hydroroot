###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to compute architecture or to construct it from text file
#   No hydraulic calculation done
#   Display root with a color according to there order
###############################################################################



import argparse
import time
import glob
import pandas as pd

from openalea.plantgl.algo.view import view
import openalea.plantgl.all as pgl
from openalea.mtg import turtle as turt

from openalea.hydroroot.init_parameter import Parameters
from openalea.hydroroot.display import get_root_visitor
from openalea.hydroroot import radius
from openalea.hydroroot.read_file import read_archi_data
from openalea.hydroroot.generator.markov import generate_g
from openalea.hydroroot.generator.measured_root import mtg_from_aqua_data

results = {}
Jv_global = 1.0

start_time = time.time()


################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################


parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("-o", "--outputfile", help="output csv file")
args = parser.parse_args()
filename = args.inputfile
output = args.outputfile
parameter.read_file(filename)

def order_scene(g1, has_radius=True, r_base=1.e-4, r_tip=5e-5, prune=None, name=None):
    """
    It is a copy of hydroroot.display.plot but the colors for the root according their order are hardcoded
    Display the architecture in plantGL Viewer with roots colors according to there order
    The radius property may be changed for display purpose.
    The MTG g1 stay unmodified

    :param g1: MTG() - the architecture to display
    :param has_radius: Boolean (True) - True use the radius property values, calculate them otehrwise according to r_base and r_tip
    :param  r_base: float (1e-4) - if has_radius is False, the radius at the base of a root whatever its order (mm)
    :param  r_tip: float (5e-5) - if has_radius is False, the radius at the tip of a root whatever its order (mm)
    :param  prune: float (None) - distance from the base of the primary after which the root is not displayed
    :param name: string (None) - if not None, the name of the saved file
    :return:
    """
    g = g1.copy() # because we may change the radius if we want
    visitor = get_root_visitor(prune=prune)

    # changing radius just for display
    r_base, r_tip = float(r_base), float(r_tip)
    if not has_radius:
        radius.discont_radius(g,r_base=r_base, r_tip=r_tip)

    turtle = turt.PglTurtle()
    turtle.down(180)
    scene = turt.TurtleFrame(g, visitor=visitor, turtle=turtle, gc=False)

    c = {}
    c[0] = [0,0,255]
    c[1] = [0,127,0]
    c[2] = [255,0,0]
    c[3] = c[2]
    for v in g.vertices_iter(scale = g.max_scale()):
        o = g.property('order')[v]
        g.property('color')[v] = c[o]

    # F. Bauget 2022-07-27: python 2 to 3
    shapes = scene.todict()
    colors = g.property('color')
    for vid in colors:
        if vid in shapes:
            for sh in shapes[vid]:
                sh.appearance =pgl.Material(colors[vid])
    scene = pgl.Scene([sh for shid in shapes.values() for sh in shid ])

    return scene
        
if __name__ == '__main__':

    filename = (glob.glob(parameter.archi['input_dir'] + parameter.archi['input_file'][0]))

    df = read_archi_data(filename[0])
    g = mtg_from_aqua_data(df, parameter.archi['segment_length'])

    # g has radius, here we set fictive radii just for visual comfort
    alpha = 0.2  # radius in millimeter identical for all orders
    s = order_scene(g, has_radius = False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4)
    view(s)

    dseeds = pd.read_csv('data_figures/sup-fig-5-B.csv')

    for id in dseeds.index:
        seed = int(dseeds.seed[id]) # otherwise this a numpy.int64 and the function random.seed does not like it
        primary_length = dseeds.primary_length[id]
        delta = dseeds.delta[id]
        nude_length = dseeds.nude_length[id]

        g =g = generate_g(seed, parameter.archi['length_data'],
                       parameter.archi['branching_variability'], delta,
                       nude_length, primary_length, parameter.archi['segment_length'],
                       parameter.archi['order_max'])
        # compute length property and parametrisation
        g = radius.compute_length(g, parameter.archi['segment_length'])

        # g has radius, here we set fictive radii just for visual comfort
        alpha = 0.2  # radius in millimeter identical for all orders
        s = order_scene(g, has_radius = False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4)
        view(s)