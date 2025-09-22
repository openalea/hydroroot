Modeling principles
===================

HydroRoot is dedicated to the modelization of hydraulic of root system architecture. It has two solvers each one having
its specificity regarding transport equations and numerical resolution:

- a water transport solver
- a solute and water transport solver

The water transport is modeling water flow with water potential difference as the only driving force. By electrical network analogy, the resolution is
done by MTG traversal and successive equivalent resistance calculation.

The solute and water transport solver considers hydrostatic and osmotic driving forces. The previous numerical resolution
is no longer valid due to coupled equations, and the system is solved in matrix form.

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   ./user/modeling_water.rst
   ./user/modeling_solutes.rst
   ./user/modeling_bc.rst

