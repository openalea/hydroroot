
# This file has been generated at Mon Feb 13 15:21:48 2023

from openalea.core import *


__name__ = 'HydroRoot'

__editable__ = True
__version__ = '1.0.0'
__license__ = 'CeCILL-C'
__authors__ = ''
__institutes__ = None
__description__ = ''
__url__ = 'http://openalea.github.io'
__icon__ = ''
__alias__ = []


__all__ = ['hydro_discont_radius', 'hydro_ordered_radius', 'hydro_linear', 'hydro_markov_binary_tree', 'hydro_shuffle_axis', 'hydro_compute_K', 'hydro_compute_k', 'hydro_fit_length', 'hydro_flux', 'hydro_plot_old', 'hydro_mtg_scene', 'hydro_colorbar', 'hydro_compute_length', 'hydro_compute_surface', 'hydro_compute_volume', 'hydro_compute_relative_position', 'hydro_fit_K', 'hydro_fit_property_from_csv', 'hydro_readCSVFile', 'hydro_read_archi_data', 'hydro_mtg_from_aqua_data', 'mtg_out_flux_mtg_out_flux']



hydro_discont_radius = Factory(name='discont_radius',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='discont_radius',
                inputs=[{'name': 'g'},
                        {'name': 'r_base', 'interface': 'IFloat', 'value': 0.0001},
                        {'name': 'r_tip', 'interface': 'IFloat', 'value': 5e-05}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_ordered_radius = Factory(name='ordered_radius',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='ordered_radius',
                inputs=[{'name': 'g'}, {'name': 'ref_radius', 'interface': 'IFloat', 'value': 0.0001}, {'name': 'order_decrease_factor', 'interface': 'IFloat', 'value': 0.75}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_linear = Factory(name='linear root',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='linear',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_markov_binary_tree = Factory(name='markov binary tree',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='markov_binary_tree',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_shuffle_axis = Factory(name='shuffle axis',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='shuffle_axis',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_compute_K = Factory(name='compute axial conductance',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='compute_K',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_compute_k = Factory(name='compute radial conductance',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='compute_k',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_fit_length = Factory(name='fit length',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='fit_length',
                inputs=[{'name': 'csvdata', 'interface': 'IFileStr'}, {'name': 'length', 'interface': 'IFloat', 'value': 0.0001}, {'name': 'k', 'interface': 'IInt', 'value': 1}, {'name': 's', 'interface': 'IFloat', 'value': 0.0}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_flux = Factory(name='flux',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='flux',
                inputs=[{'name': 'g'}, {'name': 'Jv', 'interface': 'IFloat', 'unit': 'microL/s', 'value': 0.1}, {'name': 'Pe', 'interface': 'IFloat', 'unit': 'MPa', 'value': 0.4}, {'name': 'Pbase', 'interface': 'IFloat', 'unit': 'MPa', 'value': 0.101325}, {'name': 'invert_model', 'interface': 'IBool', 'value': True}, {'name': 'k', 'interface': 'IStr', 'value': None, 'hide': True}, {'name': 'K', 'interface': 'IStr', 'value': None, 'hide': True}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_plot_old = Factory(name='plot root',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='plot_old',
                inputs=[{'name': 'g'}, {'name': 'has_radius', 'interface': 'IBool', 'value': True, 'hide': True}, {'name': 'r_base', 'interface': 'IFloat', 'value': 0.0001}, {'name': 'r_tip', 'interface': 'IFloat', 'value': 5e-05}, {'name': 'visitor', 'interface': 'IFunction'}, {'name': 'property', 'interface': 'IStr', 'value': 'radius'}, {'name': 'colormap', 'interface': 'IStr', 'value': 'jet'}, {'name': 'lognorm', 'interface': 'IBool', 'value': True}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_mtg_scene = Factory(name='root scene',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='mtg_scene',
                inputs=[{'name': 'g'}, {'name': 'has_radius', 'interface': 'IBool', 'value': True, 'hide': True}, {'name': 'r_base', 'interface': 'IFloat', 'value': 0.0001}, {'name': 'r_tip', 'interface': 'IFloat', 'value': 5e-05}, {'name': 'visitor', 'interface': 'IFunction'}, {'name': 'property', 'interface': 'IStr', 'value': 'radius'}, {'name': 'colormap', 'interface': 'IStr', 'value': 'jet'}, {'name': 'lognorm', 'interface': 'IBool', 'value': False}, {'name': 'colormap minimum', 'interface': 'IStr', 'value': None}, {'name': 'colormap maximum', 'interface': 'IStr', 'value': None}, {'name': 'prune', 'interface': 'IStr', 'value': None, 'hide': True}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_colorbar = Factory(name='colorbar',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='colorbar',
                inputs=[{'name': 'g'}, {'name': 'property_name', 'interface': 'IStr'}, {'name': 'colormap', 'interface': 'IStr', 'value': 'spectral'}, {'name': 'lognorm', 'interface': 'IBool', 'value': True}],
                outputs=[{'name': 'g'}, {'name': 'colorbar'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_compute_length = Factory(name='compute_length',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='compute_length',
                inputs=[{'name': 'g'}, {'name': 'length', 'interface': 'IStr', 'value': '1.e-4'}],
                outputs=[{'name': 'g'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_compute_surface = Factory(name='compute_surface',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='compute_surface',
                inputs=[{'name': 'g'}],
                outputs=[{'name': 'g'}, {'name': 'surface'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_compute_volume = Factory(name='compute_volume',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='compute_volume',
                inputs=[{'name': 'g'}],
                outputs=[{'name': 'g'}, {'name': 'volume'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_compute_relative_position = Factory(name='compute_relative_position',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='compute_relative_position',
                inputs=[{'name': 'g'}],
                outputs=[{'name': 'g'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_fit_K = Factory(name='fit_K',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='fit_K',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




hydro_fit_property_from_csv = Factory(name='fit_property_from_csv',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='fit_property_from_csv',
                inputs=[{'name': 'g'}, {'name': 'csvdata'}, {'name': 'prop_in', 'interface': 'IStr'}, {'name': 'prop_out', 'interface': 'IStr'}, {'name': 'smoothing degree', 'interface': 'IFloat', 'value': 1.0}, {'name': 'smoothing factor', 'interface': 'IFloat', 'value': 0.0}, {'name': 'plot', 'interface': 'IBool', 'value': True}, {'name': 'direct_input'}],
                outputs=[{'name': 'g'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_readCSVFile = Factory(name='readCSVFile',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='readCSVFile',
                inputs=[{'interface': IFileStr, 'name': 'file', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'data', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_read_archi_data = Factory(name='read_archi_data',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='read_archi_data',
                inputs=[{'interface': IFileStr, 'name': 'file', 'value': None, 'desc': ''}],
                outputs=[{'name': 'dataframe'}],
                widgetmodule=None,
                widgetclass=None,
               )




hydro_mtg_from_aqua_data = Factory(name='mtg_from_aqua_data',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='hydro',
                nodeclass='mtg_from_aqua_data',
                inputs=[{'name': 'df'}, {'name': 'vertex length', 'interface': 'IFloat', 'unit': 'm', 'value': 0.0001}],
                outputs=[{'name': 'g'}],
                widgetmodule=None,
                widgetclass=None,
               )




mtg_out_flux_mtg_out_flux = Factory(name='mtg_out_flux',
                authors=' (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='mtg_out_flux',
                nodeclass='mtg_out_flux',
                inputs=[{'name': 'g', 'interface': None, 'value': None, 'desc': ''}],
                outputs=[{'name': 'j', 'interface': IFloat, 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




