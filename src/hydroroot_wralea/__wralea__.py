
# This file has been generated at Fri Jun  1 10:57:35 2012

from openalea.core import *


__name__ = 'HydroRoot'

__editable__ = True
__description__ = ''
__license__ = 'CeCILL-C'
__url__ = 'http://openalea.gforge.inria.fr'
__alias__ = []
__version__ = '1.0.0'
__authors__ = ''
__institutes__ = None
__icon__ = ''


__all__ = []


r1=Factory(name='discont_radius',
#nodemodule='hydroroot.radius',
           nodemodule='hydro',
           nodeclass='discont_radius',
           inputs=[dict(name='g'), 
                   dict(name='r_base', interface='IStr', value='1e-4'),
                   dict(name='r_tip', interface='IStr', value='5e-5'), ],
            )
__all__.append('r1')

r2=Factory(name='ordered_radius',
#nodemodule='hydroroot.radius',
           nodemodule='hydro',
           nodeclass='ordered_radius',
           inputs=[dict(name='g'), 
                   dict(name='ref_radius', interface='IStr', value='1e-4'),
                   dict(name='order_decrease_factor', interface='IStr', value='0.75'), ],
            )
__all__.append('r2')


lr=Factory(name='linear root',
           nodemodule='hydro',
           nodeclass='linear',
            )
__all__.append('lr')

mbt=Factory(name='markov binary tree',
           nodemodule='hydro',
           nodeclass='markov_binary_tree',
            )
__all__.append('mbt')

sa=Factory(name='shuffle axis',
           nodemodule='hydro',
           nodeclass='shuffle_axis',
            )
__all__.append('sa')

cck=Factory(name='compute axial conductance',
           nodemodule='hydro',
           #nodemodule='hydroroot.conductance',
           nodeclass='compute_K',
            )
__all__.append('cck')

ccrk=Factory(name='compute radial conductance',
           nodemodule='hydro',
           #nodemodule='hydroroot.conductance',
           nodeclass='compute_k',
            )
__all__.append('ccrk')

_fitlen=Factory(name='fit length',
           nodemodule='hydro',
           nodeclass='fit_length',
            inputs=[dict(name='csvdata', interface='IFileStr'), 
                   dict(name='length', interface='IStr', value='1.e-4'),
                   dict(name='k', interface='IInt', value=1),
                   dict(name='s', interface='IFloat', value=0.),
                    ],
            )
__all__.append('_fitlen')

_flux=Factory(name='flux',
           nodemodule='hydro',
           #nodemodule='hydroroot.flux',
           nodeclass='flux',
            )
__all__.append('_flux')

_plot=Factory(name='plot root',
              nodemodule='hydro',
              nodeclass='plot_old',
           inputs=[dict(name='g'), 
                   dict(name='has_radius', interface='IBool', value=True, hide=True),
                   dict(name='r_base', interface='IStr', value='1e-4'),
                   dict(name='r_tip', interface='IStr', value='5e-5'),
                   dict(name='visitor', interface='IFunction'),
                   dict(name='property', interface='IStr', value='radius'),
                   dict(name='colormap', interface='IStr', value='spectral'),
                   dict(name='lognorm', interface='IBool', value=True),
                    ],

            )
__all__.append('_plot')

_cb=Factory(name='colorbar',
           nodemodule='hydro',
           nodeclass='colorbar',
           inputs = [dict(name='g'),
                   dict(name='property_name', interface='IStr'),
                   dict(name='colormap', interface='IStr', value='spectral'),
                   dict(name='lognorm', interface='IBool', value=True),
                   ],
            outputs = [dict(name='g'), dict(name='colorbar')],
            )
__all__.append('_cb')

cl=Factory(name='compute_length',
           nodemodule='hydro',
           nodeclass='compute_length',
           inputs = [dict(name='g'),
                   dict(name='length', interface='IStr', value='1.e-4'),
                   ],
            outputs = [dict(name='g')],
            )
__all__.append('cl')

cs=Factory(name='compute_surface',
           nodemodule='hydro',
           nodeclass='compute_surface',
           inputs = [dict(name='g'),
                   ],
            outputs = [dict(name='g'), dict(name='surface')],
            )
__all__.append('cs')

cv=Factory(name='compute_volume',
           nodemodule='hydro',
           nodeclass='compute_volume',
           inputs = [dict(name='g'),
                   ],
            outputs = [dict(name='g'), dict(name='volume')],
            )
__all__.append('cv')

crp=Factory(name='compute_relative_position',
           nodemodule='hydro',
           nodeclass='compute_relative_position',
           inputs = [dict(name='g'),
                   ],
            outputs = [dict(name='g')],
            )
__all__.append('crp')

fK=Factory(name='fit_K',
           nodemodule='hydro',
           nodeclass='fit_K',
            )
__all__.append('fK')

fpfc=Factory(name='fit_property_from_csv',
           nodemodule='hydro',
           nodeclass='fit_property_from_csv',
           inputs = [dict(name='g'),
                   dict(name='csvdata'),
                   dict(name='prop_in', interface='IStr'),
                   dict(name='prop_out', interface='IStr'),
                   dict(name='smoothing degree', value='3'),
                   dict(name='smoothing factor', value='0.'),
                   dict(name='plot', interface='IBool',value=True),
                   dict(name='direct_input'),
                   ],
            outputs = [dict(name='g')],
            )
__all__.append('fpfc')

rcf=Factory(name='readCSVFile',
           nodemodule='hydro',
           nodeclass='readCSVFile',
           #inputs=[dict(name='filename',interface: 'file')],
           # outputs=[dict(name='data')],
           inputs=[{'interface': IFileStr, 'name': 'file', 'value': None, 'desc': ''}],
           outputs=[{'interface': None, 'name': 'data', 'desc': ''}],
            )
__all__.append('rcf')

