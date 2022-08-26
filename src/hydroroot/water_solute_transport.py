import math
import numpy as np
from scipy import constants #, sparse
from scipy.sparse import csc_matrix, linalg
from openalea.mtg import traversal

def pressure_calculation(g, Temp = 298, sigma = 1.0, tau = 0.0, Ce = 0.0, Ps = 0.0, Cse = 0.0, dP = None,
                         Pe = 0.4, Pbase = 0.1, data = None, row = None, col = None, C_base = None):
    """the system of equation under matrix form is solved using a Newton-Raphson schemes, at each step a system A dx = b
    is solved by LU decomposition.

    :param g: MTG
    :param Temp: float (Default value = 298)
    :param sigma: float (Default value = 1.0)
    :param tau: float (Default value = 0.0)
    :param Ce: float (Default value = 0.0)
    :param Ps: float (Default value = 0.0)
    :param Cse: float (Default value = 0.0)
    :param dP: float (Default value = None)
    :param Pe: float (Default value = 0.4)
    :param Pbase: float (Default value = 0.1)
    :param data: numpy array (Default value = None)
    :param row: numpy array (Default value = None)
    :param col: numpy array (Default value = None)
    :param C_base: float (Default value = None)
    :returns: - g (MTG)
        - dx (array)
        - data (array)
        - row (array)
        - col (array)
    
    This function take into account the possibility to have non-permeating solute inside the root. This is, for instance,
    the case when performing cut and flow experiment in a solution containing such solute. Indeed, the non-permeating solute
    may enter the root at cut tips. It is referred to this non-permeating solute through Cpeg when C refers to the permeating solute
    penetrating radially the root.
    The presence of Cpeg may change the sap viscosity and has influence on the osmotic pressure see beginning of
    the function inside the root.

    """

    p_out = g.property('psi_out')  # it is pressure not potential
    p_in = g.property('psi_in')  # it is pressure not potential
    J_out = g.property('J_out')
    j = g.property('j')
    n = len(g) - 1
    radius = g.property('radius')
    length = g.property('length')
    b = np.zeros(3 * n)
    K = g.property('K')
    Kexp = g.property('K_exp')
    k = g.property('k')
    C = g.property('C')  # en mol/microL
    Cpeg = g.property('Cpeg')  # en mol/microL
    theta = g.property('theta')  # 1 if Pin >= Pout, 0 otherwise
    J_s = g.property('J_s')
    mu = g.property('mu')
    dKdCpeg = {}

    sigmaRT = sigma * constants.R * Temp * 1.0e3  # if C in mol/microL *1e9, and P in MPa *1e-6 => *1e3

    for v in g.vertices_iter(scale = 1):
        mu[v] = viscosity_peg(Cpeg[v], unit_factor = 8.0e6)  # mol/microL -> g/g of water
        dmu = derivative_viscosity_peg(Cpeg[v], unit_factor = 8.0e6)  # mol/microL -> g/g of water
        K[v] = 1.0 / mu[v] * Kexp[v] / length[v]
        dKdCpeg[v] = - dmu / mu[v] * K[v]

    if dP is None:
        p_ext = Pe
    else:
        p_ext = Pbase + dP
    Pi_e_peg = osmotic_p_peg(Ce, unit_factor = 8.0e6)  # from mol/microL to g/g
    psi_ext = p_ext - sigmaRT * Cse + Pi_e_peg  # subtract the osmotic pressure
    derivative = derivative_osmotic_p_peg(Ce, unit_factor = 8.0e6)  # from mol/microL to g/g

    # Select the base of the root
    v_base = 1  # it has to be = 1 here because based on first index == 1
    nid = 0
    if C_base is None: C_base = C[v_base]  # boundary condition at the root base
    m = 20 * n - 12  # -12 because of the coefficients outside the matrix at the boundaries
    ############
    # row and col indexes should be calculated once
    ############
    if (row is None) | (col is None):
        row = np.empty(m)
        col = np.empty(m)
        for v in g.vertices_iter(scale = 1):
            kids = g.children(v)
            parent = g.parent(v)

            if v != v_base:
                # dGp_i/dP_p
                row[nid] = int(3 * v - 3)
                col[nid] = int(3 * parent - 3)
                nid += 1

                # dGc_i/dP_p
                row[nid] = int(3 * v - 2)
                col[nid] = int(3 * parent - 3)
                nid += 1

                # dGCpeg_i/dP_p
                row[nid] = int(3 * v - 1)
                col[nid] = int(3 * parent - 3)
                nid += 1

                # dGc_i/dC_p
                row[nid] = int(3 * v - 2)
                col[nid] = int(3 * parent - 2)
                nid += 1

                # dGCpeg_i/dCpeg_p
                row[nid] = int(3 * v - 1)
                col[nid] = int(3 * parent - 1)
                nid += 1

            # dGp_i/dP_i
            row[nid] = int(3 * v - 3)
            col[nid] = int(3 * v - 3)
            nid += 1

            # dGc_i/dP_i
            row[nid] = int(3 * v - 2)
            col[nid] = int(3 * v - 3)
            nid += 1

            # dGCpeg_i/dP_i
            row[nid] = int(3 * v - 1)
            col[nid] = int(3 * v - 3)
            nid += 1

            # dGp_i/dC_i
            row[nid] = int(3 * v - 3)
            col[nid] = int(3 * v - 2)
            nid += 1

            # dGc_i/dC_i
            row[nid] = int(3 * v - 2)
            col[nid] = int(3 * v - 2)
            nid += 1

            # dGp_i/dCpeg_i
            row[nid] = int(3 * v - 3)
            col[nid] = int(3 * v - 1)
            nid += 1

            # dGc_i/dCpeg_i
            row[nid] = int(3 * v - 2)  # F. Bauget 2021-10-12
            col[nid] = int(3 * v - 1)
            nid += 1

            # dGCpeg_i/dCpeg_i
            row[nid] = int(3 * v - 1)
            col[nid] = int(3 * v - 1)
            nid += 1

            for cid in kids:
                # dGp_i/dP_j
                row[nid] = int(3 * v - 3)
                col[nid] = int(3 * cid - 3)
                nid += 1

                # dGp_i/dCpeg_j
                row[nid] = int(3 * v - 3)  # F. Bauget 2021-10-12
                col[nid] = int(3 * cid - 1)
                nid += 1

                # dGc_i/dP_j
                row[nid] = int(3 * v - 2)
                col[nid] = int(3 * cid - 3)
                nid += 1

                # dGCpeg_i/dP_j
                row[nid] = int(3 * v - 1)
                col[nid] = int(3 * cid - 3)
                nid += 1

                # dGc_i/dC_j
                row[nid] = int(3 * v - 2)
                col[nid] = int(3 * cid - 2)
                nid += 1

                # dGc_i/dCpeg_j
                row[nid] = int(3 * v - 2)  # F. Bauget 2021-10-12
                col[nid] = int(3 * cid - 1)
                nid += 1

                # dGCpeg_i/dCpeg_j
                row[nid] = int(3 * v - 1)
                col[nid] = int(3 * cid - 1)
                nid += 1

    ############
    # non-zero Jacobian terms
    ############
    nid = 0
    data = np.empty(m)
    for v in g.vertices_iter(scale = 1):
        kids = g.children(v)
        parent = g.parent(v)

        if v == v_base:
            p_parent = Pbase
            C_parent = C_base # C[v]  # dC/dx=0 => Cp=C[1] or C[1] = C_base
            # theta[v] = 1.0 # F. Bauget 2022-08-03
            Cpeg_parent = Cpeg[v]
        else:
            p_parent = p_in[parent]
            C_parent = C[parent]
            Cpeg_parent = Cpeg[parent]

            # dGp_i/dP_p
            data[nid] = -K[v]
            nid += 1

            # dGc_i/dP_p
            data[nid] = -K[v] * (theta[v] * C[v] + (1 - theta[v]) * C_parent)
            nid += 1

            # dGCpeg_i/dP_p
            data[nid] = -K[v] * (theta[v] * Cpeg[v] + (1 - theta[v]) * Cpeg_parent)
            nid += 1

            # dGc_i/dC_p
            data[nid] = (1 - theta[v]) * K[v] * (p_in[v] - p_parent)
            nid += 1

            # dGCpeg_i/dCpeg_p
            data[nid] = (1 - theta[v]) * K[v] * (p_in[v] - p_parent)
            nid += 1

        diag = K[v] + k[v]
        var = K[v] * (theta[v] * C[v] + (1 - theta[v]) * C_parent)
        var_peg = K[v] * (theta[v] * Cpeg[v] + (1 - theta[v]) * Cpeg_parent)
        var2 = 0.
        for cid in kids:
            diag += K[cid]
            var += K[cid] * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])
            var_peg += K[cid] * (theta[cid] * Cpeg[cid] + (1 - theta[cid]) * Cpeg[v])
            var2 += (1 - theta[cid]) * K[cid] * (p_in[cid] - p_in[v])
        if not kids:
            if g.label(v) == 'cut':
                theta_ext = int(p_in[v] <= p_ext)
                diag += K[v]
                var += K[v] * (theta_ext * Cse + (1 - theta_ext) * C[v])
                var_peg += K[v] * (theta_ext * Ce + (1 - theta_ext) * Cpeg[v])
                var2 += (1 - theta_ext) * K[v] * (p_ext - p_in[v])

        # dGp_i/dP_i
        data[nid] = diag
        nid += 1

        # dGc_i/dP_i
        data[nid] = var
        nid += 1

        # dGCpeg_i/dP_i
        data[nid] = var_peg
        nid += 1

        # dGp_i/dC_i
        data[nid] = -k[v] * sigmaRT
        nid += 1

        # dGc_i/dC_i
        data[nid] = theta[v] * K[v] * (p_in[v] - p_parent) - var2 + Ps * radius[v] * 2 * np.pi * length[v] * 1e9
        nid += 1

        # dGp_i/dCpeg_i
        derivative = derivative_osmotic_p_peg(Cpeg[v], unit_factor = 8.0e6)  # mol/microL -> g/g
        data[nid] = k[v] * derivative + dKdCpeg[v] * (p_in[v] - p_parent)  # F. Bauget 2021-10-12
        nid += 1

        # dGc_i/dCpeg_i
        data[nid] = dKdCpeg[v] * (p_in[v] - p_parent) * (
                    theta[v] * C[v] + (1 - theta[v]) * C_parent)  # F. Bauget 2021-10-12
        nid += 1

        # dGCpeg_i/dCpeg_i
        data[nid] = theta[v] * K[v] * (p_in[v] - p_parent) - var2 + \
                    dKdCpeg[v] * (p_in[v] - p_parent) * (
                                theta[v] * Cpeg[v] + (1 - theta[v]) * Cpeg_parent)  # F. Bauget 2021-10-12
        nid += 1

        for cid in kids:
            # dGp_i/dP_j
            data[nid] = -K[cid]
            nid += 1

            # dGp_i/dCpeg_j
            data[nid] = - dKdCpeg[cid] * (p_in[cid] - p_in[v])  # F. Bauget 2021-10-12
            nid += 1

            # dGc_i/dP_j
            data[nid] = -K[cid] * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])
            nid += 1

            # dGCpeg_i/dP_j
            data[nid] = -K[cid] * (theta[cid] * Cpeg[cid] + (1 - theta[cid]) * Cpeg[v])
            nid += 1

            # dGc_i/dC_j
            data[nid] = -theta[cid] * K[cid] * (p_in[cid] - p_in[v])
            nid += 1

            # dGc_i/dCpeg_j
            data[nid] = - dKdCpeg[cid] * (p_in[cid] - p_in[v]) * (
                        theta[cid] * C[cid] + (1 - theta[cid]) * C[v])  # F. Bauget 2021-10-12
            nid += 1

            # dGCpeg_i/dCpeg_j
            data[nid] = -theta[cid] * K[cid] * (p_in[cid] - p_in[v]) - \
                        dKdCpeg[cid] * (p_in[cid] - p_in[v]) * (
                                    theta[cid] * Cpeg[cid] + (1 - theta[cid]) * Cpeg[v])  # F. Bauget 2021-10-12
            nid += 1

        # -Gp
        Pi_peg = osmotic_p_peg(Cpeg[v], unit_factor = 8.0e6)  # from mol/microL to g/g
        b[3 * v - 3] = -(-K[v] * p_parent + diag * p_in[v] - k[v] * sigmaRT * C[v] - k[v] * psi_ext + k[v] * Pi_peg)
        # -Gc
        b[3 * v - 2] = -((theta[v] * C[v] + (1 - theta[v]) * C_parent) * K[v] * (p_in[v] - p_parent) -
                         J_s[v] * radius[v] * 2 * np.pi * length[v] + Ps * radius[v] * 2 * np.pi * length[v] * (
                                     C[v] - Cse) * 1e9)
        # -GCpeg
        b[3 * v - 1] = -((theta[v] * Cpeg[v] + (1 - theta[v]) * Cpeg_parent) * K[v] * (p_in[v] - p_parent))
        for cid in kids:
            b[3 * v - 3] += K[cid] * p_in[cid]  # += because it is -Gp
            b[3 * v - 2] += K[cid] * (p_in[cid] - p_in[v]) * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])  # idem
            b[3 * v - 1] += K[cid] * (p_in[cid] - p_in[v]) * (
                        theta[cid] * Cpeg[cid] + (1 - theta[cid]) * Cpeg[v])  # idem

        if not kids:
            if g.label(v) == 'cut':
                theta_ext = int(p_in[v] <= p_ext)
                b[3 * v - 3] += K[v] * p_ext  # += because it is -Gp
                b[3 * v - 2] += K[v] * (p_ext - p_in[v]) * (theta_ext * Cse + (1 - theta_ext) * C[v])  # idem
                b[3 * v - 1] += K[v] * (p_ext - p_in[v]) * (theta_ext * Ce + (1 - theta_ext) * Cpeg[v])  # idem

    A = csc_matrix((data, (row, col)), shape = (3 * n, 3 * n))

    solve = linalg.splu(A)
    dx = solve.solve(b)

    for v in g.vertices_iter(scale = 1):
        p_in[v] = p_in[v] + dx[3 * v - 3]
        C[v] = C[v] + dx[3 * v - 2]
        Cpeg[v] = Cpeg[v] + dx[3 * v - 1]
        if Cpeg[v] > Ce: Cpeg[v] = Ce
        if Cpeg[v] < 1.0e-20: Cpeg[v] = 1.0e-20
        if C[v] < 0.: C[v] = 0.

    for v in g.vertices_iter(scale = 1):
        parent = g.parent(v)
        if parent is None:
            assert v == v_base
            p_out[v] = Pbase
        else:
            p_out[v] = p_in[parent]

    for v in g.vertices_iter(scale = 1):
        Pi_peg = osmotic_p_peg(Cpeg[v], unit_factor = 8.0e6)  # from mol/microL to g/g
        J_out[v] = K[v] * (p_in[v] - p_out[v])
        j[v] = k[v] * (psi_ext - p_in[v] + sigmaRT * C[v] - Pi_peg)
        theta[v] = int(p_out[v] <= p_in[v])
        J_s[v] = tau

    return g, dx, data, row, col


def pressure_calculation_no_non_permeating_solutes(g, Temp = 298, sigma = 1.0, tau = 0.0, Ce = 0.0,
                                                   Ps = 0.0, Cse = 0.0, dP = None,
                                                   Pe = 0.4, Pbase = 0.1, data = None, row = None, col = None, C_base = None):
    """As :func:`~water_solute_transport.py.pressure_calculation` without non-permeating solutes

    :param g: MTG
    :param Temp: float (Default value = 298)
    :param sigma: float (Default value = 1.0)
    :param tau: float (Default value = 0.0)
    :param Ps: float (Default value = 0.0)
    :param Cse: float (Default value = 0.0)
    :param dP: float (Default value = None)
    :param Pe: float (Default value = 0.4)
    :param Pbase: float (Default value = 0.1)
    :param data: numpy array (Default value = None)
    :param row: numpy array (Default value = None)
    :param col: numpy array (Default value = None)
    :param C_base: float (Default value = None)
    :param Ce:  (Default value = 0.0)
    :returns: - g (MTG)
        - dx (array)
        - data (array)
        - row (array)
        - col (array)

    """

    p_out = g.property('psi_out')  # it is pressure not potential
    p_in = g.property('psi_in')  # it is pressure not potential
    J_out = g.property('J_out')
    j = g.property('j')
    n = len(g) - 1
    radius = g.property('radius')
    length = g.property('length')
    b = np.zeros(2 * n)
    K = g.property('K')
    k = g.property('k')
    C = g.property('C')  # en mol/microL
    theta = g.property('theta')  # 1 if Pin >= Pout, 0 otherwise
    J_s = g.property('J_s')

    sigmaRT = sigma * constants.R * Temp * 1.0e3  # if C in mol/microL *1e9, and P in MPa *1e-6 => *1e3

    if dP is None:
        p_ext = Pe
    else:
        p_ext = Pbase + dP
    Pi_e_peg = osmotic_p_peg(Ce, unit_factor = 8.0e6)  # from mol/microL to g/g
    psi_ext = p_ext - sigmaRT * Cse + Pi_e_peg  # subtract the osmotic pressure

    # Select the base of the root
    v_base = 1  # it has to be = 1 here because based on firts index == 1
    if C_base is None: C_base = C[v_base] # boundary condition at the root base
    nid = 0
    if data is None:
        # constant Matrix elements calculated only once, is data is None not a lot faster
        # only useful for calculation without Cpeg otherwise K changes
        m = 10 * n - 6  # -6 because of the coefficients outside the matrix at the boundaries
        data = np.empty(m)
        row = np.empty(m)
        col = np.empty(m)
        for v in g.vertices_iter(scale = 1):
            kids = g.children(v)
            parent = g.parent(v)

            if v == v_base:
                p_parent = Pbase
            else:
                p_parent = p_in[parent]

                # dGp_i/dP_p
                data[nid] = -K[v]
                row[nid] = int(2 * v - 2)
                col[nid] = int(2 * parent - 2)
                nid += 1

            diag = K[v] + k[v]
            for cid in kids:
                diag += K[cid]
            if not kids:
                if g.label(v) == 'cut':
                    diag += K[v]

            # dGp_i/dP_i
            data[nid] = diag
            row[nid] = int(2 * v - 2)
            col[nid] = int(2 * v - 2)
            nid += 1

            # dGp_i/dC_i
            data[nid] = -k[v] * sigmaRT
            row[nid] = int(2 * v - 2)
            col[nid] = int(2 * v - 1)
            nid += 1

            for cid in kids:
                # dGp_i/dP_j
                data[nid] = -K[cid]
                row[nid] = int(2 * v - 2)
                col[nid] = int(2 * cid - 2)
                nid += 1

        for v in g.vertices_iter(scale = 1):
            # indexes of Matrix elements and right hand side elements depending on C and P
            kids = g.children(v)
            parent = g.parent(v)

            if v != v_base:
                # dGc_i/dP_p
                row[nid] = int(2 * v - 1)
                col[nid] = int(2 * parent - 2)
                nid += 1

                # dGc_i/dC_p
                row[nid] = int(2 * v - 1)
                col[nid] = int(2 * parent - 1)
                nid += 1

            # dGc_i/dP_i
            row[nid] = int(2 * v - 1)
            col[nid] = int(2 * v - 2)
            nid += 1

            # dGc_i/dC_i
            row[nid] = int(2 * v - 1)
            col[nid] = int(2 * v - 1)
            nid += 1

            for cid in kids:
                # dGc_i/dP_j
                row[nid] = int(2 * v - 1)
                col[nid] = int(2 * cid - 2)
                nid += 1

                # dGc_i/dC_j
                row[nid] = int(2 * v - 1)
                col[nid] = int(2 * cid - 1)
                nid += 1

    nid = 4 * n - 2  # 2 elements outside boundaries
    for v in g.vertices_iter(scale = 1):
        # Matrix elements and right hand side elements depending on C and P
        kids = g.children(v)
        parent = g.parent(v)

        if v == v_base:
            p_parent = Pbase
            C_parent = C_base # C[v]  # dC/dx=0 => Cp=C[1] or C[1] = C_base
            # theta[v] = 1.0 # F. Bauget 2022-08-03
        else:
            p_parent = p_in[parent]
            C_parent = C[parent]
            # dGc_i/dP_p
            data[nid] = -K[v] * (theta[v] * C[v] + (1 - theta[v]) * C_parent)
            nid += 1

            # dGc_i/dC_p
            data[nid] = (1 - theta[v]) * K[v] * (p_in[v] - p_parent)
            nid += 1

        diag = K[v] + k[v]
        var = K[v] * (theta[v] * C[v] + (1 - theta[v]) * C_parent)
        var2 = 0.
        for cid in kids:
            diag += K[cid]
            var += K[cid] * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])
            var2 += (1 - theta[cid]) * K[cid] * (p_in[cid] - p_in[v])
        if not kids:
            if g.label(v) == 'cut':
                theta_ext = int(p_in[v] <= p_ext)
                diag += K[v]
                var += K[v] * (theta_ext * Cse + (1 - theta_ext) * C[v])
                var2 += (1 - theta_ext) * K[v] * (p_ext - p_in[v])

        # dGc_i/dP_i
        data[nid] = var
        nid += 1

        # dGc_i/dC_i
        data[nid] = theta[v] * K[v] * (p_in[v] - p_parent) - var2 + Ps * radius[v] * 2 * np.pi * length[v] * 1e9
        nid += 1

        for cid in kids:
            # dGc_i/dP_j
            data[nid] = -K[cid] * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])
            nid += 1

            # dGc_i/dC_j
            data[nid] = -theta[cid] * K[cid] * (p_in[cid] - p_in[v])
            nid += 1

        # -Gp
        b[2 * v - 2] = -(-K[v] * p_parent + diag * p_in[v] - k[v] * sigmaRT * C[v] - k[v] * psi_ext)
        # -Gc
        b[2 * v - 1] = -((theta[v] * C[v] + (1 - theta[v]) * C_parent) * K[v] * (p_in[v] - p_parent) -
                         J_s[v] * radius[v] * 2 * np.pi * length[v] + Ps * radius[v] * 2 * np.pi * length[v] * (
                                     C[v] - Cse) * 1e9)
        for cid in kids:
            b[2 * v - 2] += K[cid] * p_in[cid]  # += because it is -Gp
            b[2 * v - 1] += K[cid] * (p_in[cid] - p_in[v]) * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])  # idem

        if not kids:
            if g.label(v) == 'cut':
                theta_ext = int(p_in[v] <= p_ext)
                b[2 * v - 2] += K[v] * p_ext  # += because it is -Gp
                b[2 * v - 1] += K[v] * (p_ext - p_in[v]) * (theta_ext * Cse + (1 - theta_ext) * C[v])  # idem

    A = csc_matrix((data, (row, col)), shape = (2 * n, 2 * n))

    solve = linalg.splu(A)
    dx = solve.solve(b)

    for v in g.vertices_iter(scale = 1):
        p_in[v] = p_in[v] + dx[2 * v - 2]
        C[v] = C[v] + dx[2 * v - 1]
        if C[v] < 0.: C[v] = 0.

    for v in g.vertices_iter(scale = 1):
        parent = g.parent(v)
        if parent is None:
            assert v == v_base
            p_out[v] = Pbase
        else:
            p_out[v] = p_in[parent]

    for v in g.vertices_iter(scale = 1):
        J_out[v] = K[v] * (p_in[v] - p_out[v])
        theta[v] = int(p_out[v] <= p_in[v])
        J_s[v] = tau
        j[v] = k[v] * (psi_ext - p_in[v] + sigmaRT * C[v])

    return g, dx, data, row, col

def pressure_calculation_drag(g, Temp = 298, sigma = 1.0, tau=0.0, Ce = 0.0, Ps = 0.0, Cse = 0.0, dP = None,
                              Pe = 0.4, Pbase = 0.1, data = None, row = None, col = None, C_base = None):
    """Deprecated (do not work with cut and flow)
    As :func:`~water_solute_transport.py.pressure_calculation` without non-permeating solutes and with a drag term
    in the solute transport equation.

    :param g: MTG
    :param Temp: float (Default value = 298)
    :param sigma: float (Default value = 1.0)
    :param tau: float (Default value = 0.0)
    :param Ps: float (Default value = 0.0)
    :param Cse: float (Default value = 0.0)
    :param dP: float (Default value = None)
    :param Pe: float (Default value = 0.4)
    :param Pbase: float (Default value = 0.1)
    :param data: numpy array (Default value = None)
    :param row: numpy array (Default value = None)
    :param col: numpy array (Default value = None)
    :param C_base: float (Default value = None)
    :param Ce:  (Default value = 0.0)
    :returns: - g (MTG)
        - dx (array)
        - data (array)
        - row (array)
        - col (array)

    """
    # F. Bauget 2021-03-12 : added drag term -(1-sigma) * 0.5 * (C+Cse) * j
    p_out = g.property('psi_out') # it is pressure not potential
    p_in = g.property('psi_in') # it is pressure not potential
    J_out = g.property('J_out')
    j = g.property('j')
    n = len(g) - 1
    radius = g.property('radius')
    length = g.property('length')
    b = np.zeros(2*n)
    K = g.property('K')
    k = g.property('k')
    C = g.property('C') # en mol/microL
    theta = g.property('theta') # 1 if Pin >= Pout, 0 otherwise

    sigmaRT = sigma * constants.R * Temp * 1.0e3 # if C in mol/microL *1e9, and P in MPa *1e-6 => *1e3
    if dP is None:
        p_ext = Pe
    else:
        p_ext = Pbase + dP
    Pi_e_peg = osmotic_p_peg(Ce, unit_factor = 8.0e6) # from mol/microL to g/g
    psi_ext = p_ext - sigmaRT * Cse + Pi_e_peg # subtract the osmotic pressure

    # Select the base of the root
    v_base = 1 # it has to be = 1 here because based on firts index == 1
    if C_base is None: C_base = C[v_base] # boundary condition at the root base
    nid = 0
    if data is None:
        #constant Matrix elements calculated only once
        m = 10 * n - 6 #-6 because of the coefficients outside the matrix at the boundaries
        data = np.empty(m)
        row = np.empty(m)
        col = np.empty(m)
        for v in g.vertices_iter(scale = 1):
            kids = g.children(v)
            parent = g.parent(v)

            if v == v_base:
                p_parent = Pbase
            else:
                p_parent = p_in[parent]

                # dGp_i/dP_p
                data[nid] = -K[v]
                row[nid] = int(2 * v - 2)
                col[nid] = int(2 * parent - 2)
                nid += 1

            diag = K[v] + k[v]
            for cid in kids:
                diag += K[cid]

            # dGp_i/dP_i
            data[nid] = diag
            row[nid] = int(2 * v - 2)
            col[nid] = int(2 * v - 2)
            nid += 1

            # dGp_i/dC_i
            data[nid] = -k[v] * sigmaRT
            row[nid] = int(2 * v - 2)
            col[nid] = int(2 * v - 1)
            nid += 1

            for cid in kids:
                # dGp_i/dP_j
                data[nid] = -K[cid]
                row[nid] = int(2 * v - 2)
                col[nid] = int(2 * cid - 2)
                nid += 1


        for v in g.vertices_iter(scale = 1):
            # indexes of Matrix elements and right hand side elements depending on C and P
            kids = g.children(v)
            parent = g.parent(v)

            if v != v_base:
                #dGc_i/dP_p
                row[nid] = int(2*v-1)
                col[nid] = int(2*parent-2)
                nid += 1

                #dGc_i/dC_p
                row[nid] = int(2*v-1)
                col[nid] = int(2*parent-1)
                nid += 1

            #dGc_i/dP_i
            row[nid] = int(2*v-1)
            col[nid] = int(2*v-2)
            nid += 1

            #dGc_i/dC_i
            row[nid] = int(2*v-1)
            col[nid] = int(2*v-1)
            nid += 1

            for cid in kids:
                #dGc_i/dP_j
                row[nid] = int(2*v-1)
                col[nid] = int(2*cid-2)
                nid += 1

                #dGc_i/dC_j
                row[nid] = int(2*v-1)
                col[nid] = int(2*cid-1)
                nid += 1


    nid = 4 * n - 2 #2 elements outside boundaries
    for v in g.vertices_iter(scale=1):
        # Matrix elements and right hand side elements depending on C and P
        kids = g.children(v)
        parent = g.parent(v)
        dS = radius[v] * 2 * np.pi * length[v]

        if v == v_base:
            p_parent = Pbase
            C_parent = C_base # C[v]  # dC/dx=0 => Cp=C[1] or C[1] = C_base
            # theta[v] = 1.0 # F. Bauget 2022-08-03
        else:
            p_parent = p_in[parent]
            C_parent = C[parent]
            #dGc_i/dP_p
            data[nid] = -K[v] * (theta[v] * C[v] + (1 - theta[v]) * C_parent)
            nid += 1

            #dGc_i/dC_p
            data[nid] = (1 - theta[v]) * K[v] * (p_in[v] - p_parent)
            nid += 1

        diag = K[v] + k[v]
        var = K[v] * (theta[v] * C[v] + (1 - theta[v]) * C_parent)
        var2 = 0.
        for cid in kids:
            diag += K[cid]
            var += K[cid] * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])
            var2 += (1 - theta[cid]) * K[cid] * (p_in[cid] - p_in[v])

        #dGc_i/dP_i
        data[nid] = var + (1.0 - sigma) * 0.5 * (C[v] + Cse) * k[v] * dS
        nid += 1

        #dGc_i/dC_i
        data[nid] = theta[v] * K[v] * (p_in[v] - p_parent) - var2 + Ps * dS * 1e9
        data[nid] -= 0.5 * (1.0 - sigma) * dS * (j[v] + k[v] * sigmaRT * (C[v] + Cse))
        nid += 1

        for cid in kids:
            #dGc_i/dP_j
            data[nid] = -K[cid] * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v])
            nid += 1

            #dGc_i/dC_j
            data[nid] = -theta[cid] * K[cid] * (p_in[cid] - p_in[v])
            nid += 1

        #-Gp
        b[2*v-2] = -(-K[v] * p_parent + diag * p_in[v] - k[v] * sigmaRT * C[v] - k[v] * psi_ext)
        #-Gc
        b[2*v-1] = -((theta[v] * C[v] + (1 - theta[v]) * C_parent) * K[v] * (p_in[v] - p_parent) -
                     tau * dS + Ps * dS * (C[v] - Cse) * 1e9
                     - (1.0 - sigma) * 0.5 * (C[v] + Cse) * j[v] * dS)
        for cid in kids:
            b[2*v-2] += K[cid] * p_in[cid] # += because it is -Gp
            b[2*v-1] += K[cid] * (p_in[cid] - p_in[v]) * (theta[cid] * C[cid] + (1 - theta[cid]) * C[v]) #idem

    A = csc_matrix((data, (row, col)), shape=(2*n, 2*n))

    solve = linalg.splu(A)
    dx = solve.solve(b)

    for v in g.vertices_iter(scale=1):
        p_in[v] = p_in[v] + dx[2*v-2]
        C[v] = C[v] + dx[2*v-1]


    for v in g.vertices_iter(scale=1):
        parent = g.parent(v)
        if parent is None:
            assert v == v_base
            p_out[v] = Pbase
        else:
            p_out[v] = p_in[parent]

    for v in g.vertices_iter(scale=1):
        J_out[v] = K[v] * (p_in[v] - p_out[v])
        j[v] = k[v] * (psi_ext - p_in[v] + sigmaRT * C[v])
        theta[v] = int(p_out[v] <= p_in[v])

    return g, dx, data, row, col

def init_some_MTG_properties(g, tau = 0., Cini = 0., Cpeg_ini = 0., t = 1):
    """initialization of some g properties specific to solute transport:
        - 'C': the permeating solute concentration,
        - 'Cpeg': the non-permeating solute concentration,
        - 'theta':  1 or 0 depending on the local flow direction in the xylem vessels 1 from tip to base, 0 the opposite
        - 'J_s': the pumping rate

    :param g: MTG
    :param tau: float (Default value = 0.)
    :param Cini: float (Default value = 0.)
    :param Cpeg_ini: float (Default value = 0.)
    :param t: int (Default value = 1)
    :returns: - g

    """
    C = g.property('C')
    Cpeg = g.property('Cpeg')
    theta = g.property('theta') # 1 if Pin >= Pout, 0 otherwise
    J_s = g.property('J_s')

    # Select the base of the root
    v_base = next(g.component_roots_at_scale_iter(g.root, scale = g.max_scale()))
    for v in traversal.post_order2(g, v_base):
        C[v] = Cini # in mol/microL
        Cpeg[v] = Cpeg_ini
        theta[v] = t
        J_s[v] = tau

    return g

def viscosity_peg(Cpeg = 0.0, unit_factor = 1.0):
    """Dynamic viscosity, mu, calculation in mPa.s according to the PEG-8000 concentration for a temperature T = 298 K
    Fit done on a point at Cpeg=0, mu = 1 mPa.s and the 2 fist data points (Cpeg=0.1 and 0.2) of
    table 2 from Gonzllez-Tello J. Chem. Eng. Data 1994,39, 611-614
    
    mu = -17.4 + 18.4 exp(w/0.279)
    mu in mPa.s
    w in g/g water

    :param Cpeg: Float (Default value = 0.0)
    :param unit_factor: Float (Default value = 1.0)
    :returns: - mu: Float, the dynamic viscosity in mPa.s

    """

    if Cpeg <= 1.0e-20: Cpeg = 0.0
    w = Cpeg * unit_factor
    if w > 1.: w=1.0
    mu = -17.4 + 18.4 * math.exp(w/0.279) # 25 C
    # mu = -9.81 + 10.8 * math.exp(w/0.176) # 20 C

    return mu

def derivative_viscosity_peg(Cpeg = 0.0, unit_factor = 1.0):
    """Dynamic viscosity derivative according to concentration, dmu
    see viscosity_peg
    
    dmu/dw = 18.4 exp(w/0.279)/0.279
    dmu in mPa.s
    w in g/g water

    :param Cpeg: Float (Default value = 0.0)
    :param unit_factor: Float (Default value = 1.0)
    :returns: - dmu: Float, the derivative

    """

    if Cpeg <= 1.0e-20: Cpeg = 0.0
    w = Cpeg * unit_factor
    if w > 1.: w=1.0
    dmu = 18.4 * math.exp(w/0.279) / 0.279 # 25 C
    # dmu = 10.8 * math.exp(w/0.176)/0.176 # 20 C

    return dmu

def osmotic_p_peg(Cpeg = 0.0, unit_factor = 1.0, T = 25.0):
    """osmotic pressure calculation of the PEG 8000
    equation 1 Michel, Plant Physiol. (1983) 72, 66-70
    
    pi = (1.29*w^2*T - 140*w^2 - 4.0*w)
    pi: bar
    w: g/g
    T: Celcius

    :param Cpeg: Float (Default value = 0.0)
    :param unit_factor: Float (Default value = 1.0)
    :param T: Float (Default value = 25.0)
    :returns: - osmotic_p: Float, the osmotic pressure (MPa)

    """
    if Cpeg < 1.0e-20: Cpeg = 1.0e-20
    w = Cpeg * unit_factor
    osmotic_p = (1.29 * w ** 2 * T - 140 * w ** 2 - 4.0 * w) * 0.1 # bar -> MPa

    return osmotic_p

def derivative_osmotic_p_peg(Cpeg = 0.0, unit_factor = 1.0, T = 25.0):
    """the 1st derivative according to X of the osmotic pressure calculation of the PEG 8000
    equation 1 Michel, Plant Physiol. (1983) 72, 66-70
    
    dpi/dw = (2 * 1.29*w*T - 2*140*w - 4.0)
    pi: bar
    w: g/g
    T: Celcius

    :param Cpeg: Float (Default value = 0.0)
    :param unit_factor: Float (Default value = 1.0)
    :param T: Float (Default value = 25.0)
    :returns: - osmotic_p: Float, the osmotic pressure (MPa)

    """
    if Cpeg < 1.0e-20: Cpeg = 1.0e-20
    # dp/dCpeg = dw/dCpeg * dp/dw = unit_factor * dp/dw
    w = Cpeg * unit_factor # (g/g), if [Cpeg] = mol/microL then 1 g/g = 8e6 mol/microL
    osmotic_p = (2.58 * w  * T - 280 * w  - 4.0 ) * 0.1 * unit_factor # bar/(g/g) -> MPa/(Cpeg unit)

    return osmotic_p
