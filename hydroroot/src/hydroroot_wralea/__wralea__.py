
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
           nodemodule='hydro',
           nodeclass='discont_radius',
           inputs=[dict(name='g'), 
                   dict(name='r_base', interface='IFloat', value=1e-4),
                   dict(name='r_tip', interface='IFloat', value=5e-5), ],
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
           nodeclass='compute_k',
            )
__all__.append('cck')

ccrk=Factory(name='compute radial conductance',
           nodemodule='hydro',
           nodeclass='compute_K',
            )
__all__.append('ccrk')

_flux=Factory(name='flux',
           nodemodule='hydro',
           nodeclass='flux',
            )
__all__.append('_flux')

_plot=Factory(name='plot root',
              nodemodule='hydro',
              nodeclass='plot',
           inputs=[dict(name='g'), 
                   dict(name='length', interface='IFloat', value=1e-4),
                   dict(name='has_radius', interface='IBool', value=True, hide=True),
                   dict(name='r_base', interface='IFloat', value=1e-4),
                   dict(name='r_tip', interface='IFloat', value=5e-5),
                   dict(name='visitor', interface='IFunction'),
                   dict(name='protperty', interface='IStr', value='radius'),
                    ],

            )
__all__.append('_plot')

