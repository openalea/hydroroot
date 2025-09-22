import numpy as np
from openalea.hydroroot.display import mtg_scene
from openalea.plantgl.algo.view import view
from ipywidgets import interact, fixed
from oawidgets.plantgl import PlantGL

def wdplot_architecture(g_cut, cut_n_flow_length=None):

    if cut_n_flow_length is None:
        cut_n_flow_length = []

    def my_view(cut: str = 'tot', prop: str = 'j', imgsize: tuple = (800, 800), perspective: bool = True,
                zoom: float = 1,
                azimuth: float = 0, elevation: float = 0, line_width = 1.0):
        list_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K', 'axial flux': 'J_out'}

        g = g_cut[cut].copy()

        keys = list(g.property('radius').keys())
        radius = np.array(list(g.property('radius').values()))
        new_radius = radius * line_width
        g.properties()['radius'] = dict(list(zip(keys, new_radius)))

        s = mtg_scene(g, prop_cmap = list_prop[prop], has_radius = True)
        return view(scene = s, imgsize = imgsize, perspective = perspective, zoom = zoom, azimuth = azimuth, elevation = elevation)

    _list = ['tot']
    for i in cut_n_flow_length:
        _list.append(str(i))

    interact(my_view, cut = _list, prop = ['radial flux', 'xylem pressure', 'axial K', 'axial flux'],
             imgsize = fixed((800, 800)),
             perspective = False, zoom = (0.01, 1), azimuth = (-180, 180), elevation = (-90, 90), line_width = (1, 5))


def wdplot_3D(g_cut, cut_n_flow_length=None):

    if cut_n_flow_length is None:
        cut_n_flow_length = []

    def my_view(cut: str = 'tot', prop: str = 'j', line_width = 1.0):
        list_prop = {'radial flux': 'j', 'xylem pressure': 'psi_in', 'axial K': 'K', 'axial flux': 'J_out'}

        g = g_cut[cut].copy()

        keys = list(g.property('radius').keys())
        radius = np.array(list(g.property('radius').values()))
        new_radius = radius * line_width
        g.properties()['radius'] = dict(list(zip(keys, new_radius)))

        s = mtg_scene(g, prop_cmap = list_prop[prop], has_radius = True)

        return PlantGL(s)

    _list = ['tot']
    for i in cut_n_flow_length:
        _list.append(str(i))

    interact(my_view, cut = _list, prop = ['radial flux', 'xylem pressure', 'axial K', 'axial flux'],
             line_width = (1, 5))