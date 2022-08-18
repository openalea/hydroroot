Implementation
--------------

The purpose is then to solve the mass balance equation for both fluids,
i.e. to solve on each grid block i:

.. math::

   G_{i} = 0\ \text{with}\ G_{i} = \begin{pmatrix}
   G_{{nr}_{i}} \\
   G_{r_{i}} \\
   \end{pmatrix}

or considering the whole grid system:

.. math::

   G = 0\ \text{with}\ G = \begin{pmatrix}
   {G_{nr}}_{1}{G_{r}}_{1} & \cdots & {G_{nr}}_{N}{G_{r}}_{N} \\
   \end{pmatrix}

The dimension of G is 2N, N grid blocks for G\ :sub:`nr` and
G\ :sub:`r`.

The system may be expressed according to two unknown the saturation of
one fluid and the pressure of one fluid. The saturation will always be
the reference fluid saturation. The pressure will depend on the
experiment type.

The unknowns are stored in a 2N vector X, N element for S and N element
for P. In order to have a band diagonal system, S and P elements are
store as follow:

.. math::

   X = \begin{pmatrix}
   S_{1}P_{1} & \cdots & S_{i}P_{i} & \cdots & S_{N}P_{N} \\
   \end{pmatrix}

Now to solve the system a Newton-Raphson is used leading to:

.. math:: J\ dX = - G

J is the Jacobian of G according to S and P, dX is the 2N vector
containing alternatively dS and dP. The full system is band diagonal
with for a given grid block i:

.. math::

   \begin{pmatrix}
   \frac{d{G_{nr}}_{i}}{dS_{i - 1}} & \frac{d{G_{nr}}_{i}}{dP_{i - 1}} & \frac{d{G_{nr}}_{i}}{dS_{i}} & \frac{d{G_{nr}}_{i}}{dP_{i}} & \frac{d{G_{nr}}_{i}}{dS_{i + 1}} & \frac{d{G_{nr}}_{i}}{dP_{i + 1}} \\
   \frac{d{G_{r}}_{i}}{dS_{i - 1}} & \frac{d{G_{r}}_{i}}{dP_{i - 1}} & \frac{d{G_{r}}_{i}}{dS_{i}} & \frac{d{G_{r}}_{i}}{dP_{i}} & \frac{d{G_{r}}_{i}}{dS_{i + 1}} & \frac{d{G_{r}}_{i}}{dP_{i + 1}} \\
   \end{pmatrix}\begin{pmatrix}
   \begin{matrix}
   \binom{dS_{i - 1}}{dP_{i - 1}} \\
   \binom{dS_{i}}{dP_{i}} \\
   \end{matrix} \\
   \binom{dS_{i + 1}}{dP_{i + 1}} \\
   \end{pmatrix} = - \begin{pmatrix}
   G_{{nr}_{i}} \\
   G_{r_{i}} \\
   \end{pmatrix}

The diagonal is formed by :math:`\frac{d{G_{nr}}_{i}}{dS_{i}}` and
:math:`\frac{d{G_{r}}_{i}}{dP_{i}}` so the maximum number of non-zero
element from the diagonal is 3, 3 on the left in the first line, 3 on
the right in the 2\ :sup:`nd` line. Therefore the system is 7 band
diagonal.

All other elements are strictly zero, so to avoid stocking null element
and performing loop over them J is stock in a compact way. Only the
seven elements of the diagonal are stocked giving an array of 2N lines
and 7 columns. This gives for the grid block i:

.. math::

   \begin{pmatrix}
   0 & \frac{d{G_{nr}}_{i}}{dS_{i - 1}} & \frac{d{G_{nr}}_{i}}{dP_{i - 1}} & \frac{d{G_{nr}}_{i}}{dS_{i}} & \frac{d{G_{nr}}_{i}}{dP_{i}} & \frac{d{G_{nr}}_{i}}{dS_{i + 1}} & \frac{d{G_{nr}}_{i}}{dP_{i + 1}} \\
   \frac{d{G_{r}}_{i}}{dS_{i - 1}} & \frac{d{G_{r}}_{i}}{dP_{i - 1}} & \frac{d{G_{r}}_{i}}{dS_{i}} & \frac{d{G_{r}}_{i}}{dP_{i}} & \frac{d{G_{r}}_{i}}{dS_{i + 1}} & \frac{d{G_{r}}_{i}}{dP_{i + 1}} & 0 \\
   \end{pmatrix}

**Remark (09/19/2012)**: For the two first lines, the two first element
of each do not need to be calculated because they are outside the
system, ditto for the two last lines and their two last elements. These
elements are actually calculated in the routines
“FullyImplicit_Calcul_DP_DS” and “FullyImplicit_Compr_Calcul_DP_DS”,
which is useless and so not efficient but harmless. To be corrected.

Now to solve the linear system JdX=-G, an LU decomposition is used see §
“LU decomposition J=LU” and § “Solving the system
:math:`\underline{\underline{A}}.\ \underline{x} = \underline{b}`\ ” in
Appendices.
