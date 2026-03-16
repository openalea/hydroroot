Water transport model
---------------------

The following lines are taken, or inspired, from Boursiac *et al.* 2022 [boursiac2022]_.

The hydraulic aspects of HydroRoot consists in two main components: the radial water flow between the external
medium and the xylem vessels and the axial transport through the xylem vessels. Following Doussan and 
colleagues (Doussan et al., 1998a; Doussan et al., 1998b), the root is discretized as a network of elementary
segments consisting of a microcircuit containing both radial (:math:`k`) and axial (:math:`K`)
hydraulic conductances. The elementary segments are represented by the MTG vertices. For instance, at the ith vertex, the
local radial flux is written as follow

.. math:: j_{i} = k_{i}\left( P_{e_{i}} - P_{i} \right)S_{i}

and the local axial flow is given by

.. math:: J_{i} = K_{i}\ \frac{\left( P_{out_{i}} - P_{i} \right)}{l_{i}}

:math:`P_{e_{i}}, P_{i}\ \text{and}\ P_{out_{i}}` are the water hydrostatic pressure of the external medium at vertex *i*, in
the xylem vessels of vertex *i* and in the xylem vessels of the vertex downstream to vertex *i*, respectively. 
:math:`S_{i}\ \text{and}\ l_{i}` are, respectively, the surface area and the length of the elementary segments. 

By analogy with Ohm’s law, the hydraulic architecture may be assimilated to an electrical network (Doussan et al., 1998a;
Prusinkiewicz et al., 2007). According to the boundary conditions (uniform pressure around the root
and atmospheric pressure at its base), we are able to calculate the equivalent resistance of the network and then calculate
the outflow rate. In brief, let us consider an elementary segment *i*, with :math:`R_{i} = L_{i}/K_{i}` and
:math:`r_{i} = 1/\left( k_{i}S_{i} \right)` as axial and radial resistances, respectively. Its equivalent resistance
:math:`{R_{eq}}_{i}` is calculated as follows, assuming that the apical equivalent resistance :math:`R_{eq_{i - 1}\ }` is known:

.. math:: {R_{eq}}_{i} = \frac{R_{eq_{i - 1}\ } r_{i}}{R_{eq_{i - 1}\ } + r_{i}} + R_{i}

By implementing this equation, step by step from the tips, and by
considering a branched root as a parallel network, we end up with an
equivalent resistance for the whole network, and as a consequence, an
equivalent hydraulic conductance :math:`K_{eq}` (Albasha et al., 2019;
Prusinkiewicz et al., 2007). The basal outgoing flux (*J*\ :sub:`v`) is
then calculated according to:

.. math:: J_{v} = K_{eq}\text{(}P_{e} - P_{base}\text{)}
