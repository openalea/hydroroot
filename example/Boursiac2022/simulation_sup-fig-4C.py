###############################################################################
# Date: 2022-08-11
# F. Bauget
#   Use of HydroRoot to compute architecture or to construct it from text file
#   perform simulation with s radial conductivity as a step function see
#   hydroroot.conductance.radial_step
#   display a 3D representation with the radial flux j
###############################################################################

import glob
import argparse
import tempfile, os

import openalea.plantgl.all as pgl
from IPython.display import Image, display

from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters
from hydroroot.display import plot
from hydroroot import radius
from hydroroot.conductance import axial, radial_step
from hydroroot.read_file import read_archi_data
from hydroroot.generator.measured_root import mtg_from_aqua_data


parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
args = parser.parse_args()
filename = args.inputfile
parameter.read_file(filename)

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    kexp_axial_data = axial(axial_data, axfold)
    k_radial_data = radial_step(k_radial, factor =  k_factor, x_step =  x_k_step, dx = parameter.archi['segment_length'], scale = radfold)

    # compute local jv and psi, global Jv, Keq
    g, keq, jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                   psi_e = parameter.exp['psi_e'],
                                   psi_base = parameter.exp['psi_base'],
                                   axial_conductivity_data = kexp_axial_data,
                                   radial_conductivity_data = k_radial_data)

    return g, keq, jv

if __name__ == '__main__':

    k_factor = 3.0 # kstep used in the paper k_factor = 3
    x_k_step = 0.02

    # files names of the architecture if reconstructed from a file
    # if not we just give a dummy name for the loop used to launch run
    filename = []

    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))

    for f in filename:
        df = read_archi_data(f)
        g = mtg_from_aqua_data(df, parameter.archi['segment_length'])

        # compute radius property on MTG
        g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])

        # compute length property and parametrisation
        g = radius.compute_length(g, parameter.archi['segment_length'])
        g = radius.compute_relative_position(g)

        for axfold in parameter.output['axfold']:
            g, Keq, Jv = hydro_calculation(g, axfold = axfold)

            index = f.replace(glob.glob(parameter.archi['input_dir'])[0],"")

            print((index, axfold))

            # g has radius, here we set fictive radii just for visual comfort
            alpha = 0.2 # radius in millimeter identical for all orders
            gcopy = g.copy()  # copy because we change the radius property in plot below
            plot(gcopy, has_radius=False, r_base = alpha * 1.e-3, r_tip = alpha * 9.9e-4, prop_cmap = 'j')

            # to display in the notebook, comment to display in the 3D viewer
            pgl.Viewer.widgetGeometry.setSize(450, 600) # set the picture size in px
            fn = tempfile.mktemp(suffix='.png')
            pgl.Viewer.saveSnapshot(fn)
            pgl.Viewer.stop()
            img = Image(fn)
            os.unlink(fn)
            display(img)
