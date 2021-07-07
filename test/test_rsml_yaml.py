
import rsml
from openalea.mtg import MTG, traversal

from hydroroot.hydro_io import export_mtg_to_rsml, import_rsml_to_discrete_mtg
from hydroroot import radius
from hydroroot.init_parameter import Parameters

def closed(delta, eps=1e-14, txt=''):
    if not txt:
        assert abs(delta) < eps
    else:
        assert abs(delta) < eps, txt

def root_creation(g, segment_length, ref_radius, order_decrease_factor):
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

    g = radius.ordered_radius(g, ref_radius, order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # Calculation of the distance from base of each vertex, used for cut and flow
    # Remark: this calculation is done in flux.segments_at_length; analysis.nb_roots but there is a concern with the
    # parameter dl which should be equal to vertex length but which is not pass
    _mylength = {}
    for v in traversal.pre_order2(g, 1):
        pid = g.parent(v)
        _mylength[v] = _mylength[pid] + segment_length if pid else segment_length
    g.properties()['mylength'] = _mylength

    # _length is the total length of the RSA (sum of the length of all the segments)
    _length = g.nb_vertices(scale = 1) * segment_length
    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    v_base = g.component_roots_at_scale_iter(g.root, scale = g.max_scale()).next()
    primary_length = g.property('position')[v_base]

    return g, primary_length, _length, surface

def test_rsml():
    # F. Bauget 2021-06-11 : test rsml import-export
    # 1st - import rsml to discrete MTG usable by HydroRoot
    # 2d - export this MTG to an rsml file, and reimport it
    # test: difference between primary length, total length and surface

    segment_length = 1.0e-4 #0.1 mm
    order_decrease_factor = 0.7
    ref_radius = 7e-05

    filename = 'data/arabidopsis-simple.rsml'

    rsml_units_to_metre = {}
    rsml_units_to_metre['m'] = 1.0
    rsml_units_to_metre['cm'] = 1.0e-2
    rsml_units_to_metre['mm'] = 1.0e-3
    rsml_units_to_metre['um'] = 1.0e-6
    rsml_units_to_metre['nm'] = 1.0e-9
    
    g_c = rsml.rsml2mtg(filename)
    
    resolution = g_c.graph_properties()['metadata']['resolution']
    unit = g_c.graph_properties()['metadata']['unit']
    resolution *= rsml_units_to_metre[unit]  # rsml file unit to meter
    g = import_rsml_to_discrete_mtg(g_c, segment_length = segment_length, resolution = resolution)
    
    # calculation of g properties: radius, mylength, etc.
    g, primary_length, _length, surface = root_creation(g, segment_length, ref_radius, order_decrease_factor)

    export_mtg_to_rsml(g, "test_rsml_io.rsml", segment_length = segment_length)
    g_c = rsml.rsml2mtg("test_rsml_io.rsml")
    resolution = g_c.graph_properties()['metadata']['resolution']
    unit = g_c.graph_properties()['metadata']['unit']
    resolution *= rsml_units_to_metre[unit] # rsml file unit to meter
    g2 = import_rsml_to_discrete_mtg(g_c, segment_length = segment_length, resolution = resolution)
    g2, primary_length2, _length2, surface2 = root_creation(g2, segment_length, ref_radius, order_decrease_factor)

    diff = abs(_length2 - _length + surface - surface2 + primary_length - primary_length2)

    closed(diff, txt = 'Exported rsml from discrete MTG does not give same length, surface and volume than the original imported rsml.')

def test_yaml():
    # F. Bauget 2021-06-11: test yaml file reading

    parameter = Parameters()
    parameter.read_file('parameters_test_yaml.yml')

    boolflag = (parameter.archi['branching_delay'] == [0.001])
    boolflag = (parameter.archi['branching_variability'] == 0.25)
    boolflag = (parameter.archi['input_dir'] == 'data/')
    boolflag = (parameter.archi['input_file'] == ['hydroroot-1ter.txt'])
    boolflag = (parameter.archi['length_file'] == ['data/length*order1*.csv', 'data/length*order2*.csv'])
    boolflag = (parameter.archi['nude_length'] == [0.005])
    boolflag = (parameter.archi['order_decrease_factor'] == 0.7)
    boolflag = (parameter.archi['order_max'] == 4)
    boolflag = (parameter.archi['primary_length'] == [13.0])
    boolflag = (parameter.archi['read_architecture'] == True)
    boolflag = (parameter.archi['ref_radius'] == 7e-05)
    boolflag = (parameter.archi['seed'] == [6973738])
    boolflag = (parameter.archi['segment_length'] == 0.0001)
    boolflag = (parameter.hydro['k0'] == 92.0)
    boolflag = (parameter.hydro['axial_conductance_data'] == [[0, 0.015, 0.04, 0.06, 0.0765, 0.0965, 0.1175, 0.1485, 0.1975], 
                                               [2.9e-05, 0.000124, 0.00177, 0.00218, 0.00183, 0.00305, 0.00419, 0.00345, 0.0026]])
    boolflag = (parameter.exp['Jv'] == 0.001)
    boolflag = (parameter.exp['psi_base'] == 0.101325)
    boolflag = (parameter.exp['psi_e'] == 0.401325)
    boolflag = (parameter.output['run_nb'] == 1)
    boolflag = (parameter.output['intercepts'] == [])
    boolflag = (parameter.output['axfold'] == [1.0])
    boolflag = (parameter.output['radfold'] == [0.5, 1.0, 1.5, 2.0])
    
    closed(abs(1.0 - float(boolflag)), txt = 'Error during yaml file reading.')