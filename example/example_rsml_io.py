###############################################################################
#
# Authors: F. Bauget
# Date : 2021-06-10
#
# Test import-export rsml
###############################################################################

import glob
import rsml
import argparse
import sys

from openalea.mtg import traversal

from hydroroot import radius
from hydroroot.main import hydroroot_flow
from hydroroot.init_parameter import Parameters  # import work in progress for reading init file
from hydroroot.hydro_io import export_mtg_to_rsml, import_rsml_to_discrete_mtg

from shared_functions import radial, axial

parameter = Parameters()

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="yaml input file")
args = parser.parse_args()
filename = args.inputfile
parameter.read_file(filename)

def root_builder(g):
    """
    Set MTG properties and perform some gemetrical calculation

    The vertex radius properties is set.
    The following properties are computed: length, position, mylength, surface, volume, total length,
        primary root length

    :param:
        - `g` (MTG)

    :return:
        `g`: MTG with the different properties set or computed (see comments above),
        `primary_length`: primary root length (m)
        `_length`: total root length (m)
        `surface`: total root surface (m^2)
    """

    g = radius.ordered_radius(g, parameter.archi['ref_radius'], parameter.archi['order_decrease_factor'])

    # compute length property and parametrisation
    g = radius.compute_length(g, parameter.archi['segment_length'])
    g = radius.compute_relative_position(g)

    # Calculation of the distance from base of each vertex, used for cut and flow
    # Remark: this calculation is done in flux.segments_at_length; analysis.nb_roots but there is a concern with the
    # parameter dl which should be equal to vertex length but which is not pass
    _mylength = {}
    for v in traversal.pre_order2(g, 1):
        pid = g.parent(v)
        _mylength[v] = _mylength[pid] + parameter.archi['segment_length'] if pid else parameter.archi['segment_length']
    g.properties()['mylength'] = _mylength

    # _length is the total length of the RSA (sum of the length of all the segments)
    _length = g.nb_vertices(scale = 1) * parameter.archi['segment_length']
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)


    v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))
    primary_length = g.property('position')[v_base]

    return g, primary_length, _length, surface

def hydro_calculation(g, axfold = 1., radfold = 1., axial_data = None, k_radial = None):
    if axial_data is None: axial_data = parameter.hydro['axial_conductance_data']
    if k_radial is None: k_radial = parameter.hydro['k0']
    # compute axial & radial
    Kexp_axial_data = axial(axial_data, axfold)
    k_radial_data = radial(k_radial, axial_data, radfold)

    # compute local jv and psi, global Jv, Keq
    g, Keq, Jv = hydroroot_flow(g, segment_length = parameter.archi['segment_length'],
                                   k0 = k_radial,
                                   Jv = parameter.exp['Jv'],
                                   psi_e = parameter.exp['psi_e'],
                                   psi_base = parameter.exp['psi_base'],
                                   axial_conductivity_data = Kexp_axial_data,
                                   radial_conductivity_data = k_radial_data)

    return g, Keq, Jv

if __name__ == '__main__':
    """"
    Script for testing rsml i/o
    run %run example_rsml_io parameters_test_rsml_io.yml
    
    - import 'data/arabidopsis-simple.rsml' file from https://rootsystemml.github.io/examples/arabidopsis_simple
      into continuous mtg representation
    - set the resolution according to the unit precised in the rsml, in Hydroroot lengths are in meter
    - transform the continuous mtg to discrete mtg usable in Hydroroot according to resolution
    - calculation of g properties (radius, mylength) and  flux
    
    - export discrete mtg to rsml
    - re-import rsml to continuous mtg, transform the latter to discrete mtg
    - re-do calculation of g properties (radius, mylength) and  flux
    
    - compare both calculations
    """
    axfold = parameter.output['axfold'][0]
    radfold = parameter.output['radfold'][0]

    rsml_units_to_metre = {}
    rsml_units_to_metre['m'] = 1.0
    rsml_units_to_metre['cm'] = 1.0e-2
    rsml_units_to_metre['mm'] = 1.0e-3
    rsml_units_to_metre['um'] = 1.0e-6
    rsml_units_to_metre['nm'] = 1.0e-9

    filename = []
    for f in parameter.archi['input_file']:
        filename = filename + (glob.glob(parameter.archi['input_dir'] + f))

    # import rsml into continuous mtg representation
    g_c = rsml.rsml2mtg(filename[0])
    resolution = g_c.graph_properties()['metadata']['resolution']
    unit = g_c.graph_properties()['metadata']['unit']

    if unit not in list(rsml_units_to_metre.keys()):
        sys.exit('wrong unit in rsml file, unit must be one of the following: m, cm, mm, um, nm.')

    resolution *= rsml_units_to_metre[unit] # rsml file unit to meter

    # continuous mtg to discrete mtg
    g = import_rsml_to_discrete_mtg(g_c, segment_length = parameter.archi['segment_length'], resolution = resolution)

    # calculation of g properties: radius, mylength, etc.
    g, primary_length, _length, surface = root_builder(g)

    # flux calculation
    g, Keq, Jv = hydro_calculation(g, axfold = axfold, radfold = radfold)

    print('water flux from rsml file is ', Jv, ' uL/s')

    # export g to rsml
    export_mtg_to_rsml(g, "example_rsml_export.rsml", segment_length = parameter.archi['segment_length'])

    # redo the import from rsml to discrete MTG from the exported one
    g_c = rsml.rsml2mtg("example_rsml_export.rsml")
    resolution = g_c.graph_properties()['metadata']['resolution']
    unit = g_c.graph_properties()['metadata']['unit']
    resolution *= rsml_units_to_metre[unit] # rsml file unit to meter
    g2 = import_rsml_to_discrete_mtg(g_c, segment_length = parameter.archi['segment_length'], resolution = resolution)
    g2, primary_length2, _length2, surface2 = root_builder(g2)
    g2, Keq2, Jv2 = hydro_calculation(g2, axfold = axfold, radfold = radfold)

    print('water flux from exported-imported to rsml MTG is ', Jv, ' uL/s')
    #
    print('difference in: primary_length, _length, surface, Keq, Jv are:', primary_length-primary_length2, _length-_length2, surface-surface2, Keq-Keq2, Jv-Jv2)

    # use plot() from shared_functions to display g and g2
