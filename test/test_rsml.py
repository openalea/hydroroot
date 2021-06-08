from openalea.mtg import *
from hydroroot import display
from path import Path

data = Path('./data/test_rsml.mtg')
g = MTG(data)

# 2 Compute 3D polylines
visitor = display.get_root_visitor_with_point()
scene = display.plot(g, visitor=visitor)

# Position is position3D

# 1 Scale insertion
def quotient_axis(v):
    if g.edge_type(v) == '+':
        return True
    elif g.parent(v) is None:
        return True
    else:
        return False

g = g.insert_scale(inf_scale=1, partition=quotient_axis, default_label='A')

root_id = next(g.component_roots_iter(g.root))
g = g.insert_scale(inf_scale=1, partition=lambda v: v==root_id, default_label='P')

#3 Transform to rsml 
from rsml import continuous, io

g = continuous.discrete_to_continuous(g, position='position3D')
io.mtg2rsml(g, 'myfile.rsml')
