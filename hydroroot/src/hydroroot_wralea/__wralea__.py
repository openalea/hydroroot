
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

_flux=Factory(name='flux',
           nodemodule='hydro',
           #nodemodule='hydroroot.flux',
           nodeclass='flux',
            )
__all__.append('_flux')

_plot=Factory(name='plot root',
              nodemodule='hydro',
              nodeclass='plot',
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
                   dict(name='prop_in', interface='IStr', value='None'),
                   dict(name='prop_out', interface='IStr', value='None'),
                   dict(name='smoothing degree', value='3'),
                   dict(name='smoothing factor', value='0.'),
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

