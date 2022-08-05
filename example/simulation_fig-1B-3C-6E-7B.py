###############################################################################
# Date: 2021-06-18
# F. Bauget
#   Use of HydroRoot to compute architecture or to construct it from texte file
#   Display them with a colormap based on a MTG property ('j', 'order', etc.)
#       argument --prop passed through command line
#   If order chosen then display root with a color according to there order
###############################################################################

######
# Imports

# VERSION = 2
# F. Bauget 2021-12-14: removed unused import when migration to python 3 was done
import glob
import argparse
import tempfile, os

import openalea.plantgl.all as pgl
from openalea.mtg import turtle as turt
from IPython.display import Image, display

from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters
from hydroroot.display import get_root_visitor, mtg_scene, plot
from hydroroot import radius
from hydroroot.generator.markov import my_seed,generate_g
from hydroroot.generator.measured_root import mtg_from_aqua_data
from hydroroot.read_file import read_archi_data
from hydroroot.conductance import axial, radial


################################################
# get the model parameters, the length laws are
# calculated from the files given in the yaml file
###############################################

parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
parser.add_argument("--prop", help="property to display, e.g.: order, j")
args = parser.parse_args()
filename = args.inputfile
prop = args.prop
if prop is None: prop = 'order'
parameter.read_file(filename)

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    kexp_axial_data = axial(axial_data, axfold)
    k_radial_data = radial(k_radial, axial_data, radfold)

    # compute local jv and psi, global Jv, Keq
    g, keq, jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                   psi_e = parameter.exp['psi_e'],
                                   psi_base = parameter.exp['psi_base'],
                                   axial_conductivity_data = kexp_axial_data,
                                   radial_conductivity_data = k_radial_data)

    return g, keq, jv

def plot_order(g1, has_radius=True, r_base=1.e-4, r_tip=5e-5, prune=None, name=None):
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

# F. Bauget 2022-03-15: WIP python 2  to 3
    shapes = scene.todict()
    colors = g.property('color')
    for vid in colors:
        if vid in shapes:
            for sh in shapes[vid]:
                sh.appearance =pgl.Material(colors[vid])
    scene = pgl.Scene([sh for shid in shapes.values() for sh in shid ])

    pgl.Viewer.display(scene)
    if name is not None:
            pgl.Viewer.frameGL.saveImage(name)

if __name__ == '__main__':


    # files names of the architecture if reconstructed from a file
    # if not we just give a dummy name for the loop used to launch run
    filename = []
    if parameter.archi['read_architecture']:
        run_nb = 1
        parameter.archi['seed'] = [1] # to give something to the For-Loop below
        for f in parameter.archi['input_file']:
            filename = filename + (glob.glob(parameter.archi['input_dir'] + f))
    else:
        run_nb = parameter.output['run_nb']
        filename = ['one_run']  # just to have one run in the For-Loop below
        if parameter.archi['seed'] is None:
            s = my_seed()
            parameter.archi['seed'] = list(s)


    for seed in parameter.archi['seed']:
        for f in filename:
            if parameter.archi['read_architecture']:
                df = read_archi_data(f)
                g = mtg_from_aqua_data(df, parameter.archi['segment_length'])
            else:
                length_data = parameter.archi['length_data']
                g = generate_g(seed, length_data,
                               parameter.archi['branching_variability'], parameter.archi['branching_delay'][0],
                               parameter.archi['nude_length'][0], parameter.archi['primary_length'][0],
                               parameter.archi['segment_length'],parameter.archi['order_max'])

            # compute radius property on MTG
            g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])

            # compute length property and parametrisation
            g = radius.compute_length(g, parameter.archi['segment_length'])
            g = radius.compute_relative_position(g)

            for axfold in parameter.output['axfold']:
                g, Keq, Jv = hydro_calculation(g, axfold = axfold)

                if parameter.archi['read_architecture']:
                    index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")
                else:
                    index = seed

                #prop_cmap: property to plot e.g.: for figure 1B prop_'order', for figure 3CD 'j'
                # it could also be 'J_out' the axial flux
                print(index, axfold)
                # prop_cmap = 'order'

                # g has radius, here we set fictive radii just for visual comfort
                alpha = 0.2 # radius in millimeter identical for all orders
                gcopy = g.copy() # copy because we change the radius property in plot below
                if prop != 'order':
                    plot(gcopy, has_radius=False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4, prop_cmap = prop)
                else:
                    plot_order(gcopy, has_radius=False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4)

                # for display in the notebook, comment to display in the 3D viewer
                pgl.Viewer.widgetGeometry.setSize(450, 600) # set the picture size in px
                fn = tempfile.mktemp(suffix='.png')
                pgl.Viewer.saveSnapshot(fn)
                pgl.Viewer.stop()
                img = Image(fn)
                os.unlink(fn)
                display(img)
