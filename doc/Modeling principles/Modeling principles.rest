.. vim: syntax=rst

Modeling principles
-------------------

The HydroRoot model was developed in a Python programming language, as a
component of the OpenAlea platform (Pradal et al., 2008; Pradal et al.,
2015). HydroRoot uses a Multiscale Tree Graph (MTG) (Godin and Caraglio,
1998) to represent root hydraulic architecture, which consists of the
topology of a root system (branching positions, root lengths, root
radii, etc.) and its hydraulic structure (local radial and axial
conductivities). The RSML format (Lobet et al., 2015) is used to import
and export the data to/from the HydroRoot model. The model is open
source and available through its public repository
(https://github.com/openalea/hydroroot).

The hydraulic aspects of HydroRoot consisted in two main components: the
radial water flow between the bathing solution and the xylem vessels and
the axial transport through the xylem vessels. Following Doussan and
colleagues (Doussan et al., 1998a; Doussan et al., 1998b), the root was
discretized as a network of elementary segments consisting of a
microcircuit containing both radial (*k\ i*) and axial (*K\ i*)
hydraulic conductances (Figure 1C). The local radial flux was written as
:math:`j_{i} = k_{i}\left( \psi_{e_{i}} - \psi_{i} \right)S_{i}\ `\ and
the local axial flow as
:math:`J_{i} = K_{i}\ \frac{\left( \psi_{out} - \psi_{i} \right)}{L_{i}}`,
S\ :sub:`i` and L\ :sub:`i` are respectively the surface area and the
length of the elementary segments. By analogy with Ohm’s law, both
1/(*k\ i\ S\ i)* and *L\ i*/*K\ i* may be modeled as electric
resistances, and the hydraulic architecture may be assimilated to an
electrical network (Doussan et al., 1998a; Prusinkiewicz et al., 2007).
According to the boundary conditions (uniform pressure around the root
and atmospheric pressure at its base), we are able to calculate the
equivalent resistance of the network and then calculate the outflow
rate. In brief, let us consider an elementary segment *i*, with
:math:`R_{i} = L_{i}/K_{i}` and
:math:`r_{i} = 1/\left( k_{i}S_{i} \right)` as axial and radial
resistances, respectively. Its equivalent resistance
:math:`{R_{eq}}_{i}` is calculated as follows, assuming that the apical
equivalent resistance :math:`R_{eq_{i - 1}\ }` is known:

.. math:: \frac{1}{{R_{eq}}_{i}} = \frac{1}{R_{eq_{i - 1}\ } + r_{i}} + \frac{1}{R_{i}}

By implementing this equation, step by step from the tips, and by
considering a branched root as a parallel network, we end up with an
equivalent resistance for the whole network, and as a consequence, an
equivalent hydraulic conductance :math:`K_{eq}` (Albasha et al., 2019;
Prusinkiewicz et al., 2007). The basal outgoing flux (*J*\ :sub:`v`) is
then calculated according to:

.. math:: J_{v} = K_{eq}\text{(}\psi_{e} - \psi_{base}\text{)}
