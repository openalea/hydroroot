"""F. Bauget 2020-08-13 :
    just some routine with data
"""


def flux_intercepts_for_optim(filename = None):
    """
    For cut and flow experiment.
    Return the flux of the full root, the fluxes and primary root length for the cuts for a given plant
    The plant name shoud be in the architecture file name.
    :param filename: the file name of the root architecture
    :return: J[plant] float: flux of the full root
            _Jv[plant] float: fluxes of the different root cuts
            cut_length[plant] float: primary root length at the different cuts
    """
    J = {}
    _Jv = {}
    cut_length = {}

    ### Arabidopsis
    J['200219-YBFB-esb11'] = 0.00284
    J['160316#2'] = 3.49E-02  # 300 kPa
    J['2020-01-19-9h06'] = 0.00456
    J['15012020-1045'] = 0.0062
    J['1045-sim'] = 0.00572270032329
    J['200219-YBFB-anac038'] = 0.0122
    J['200219-YBFB-Col'] = 0.00524  # col + offset
    J['200219-YBFB-horst1'] = 0.0116
    J['200619-YBFB-Col-1'] = 0.0163
    J['200619-YBFB-Col-2'] = 0.0342
    J['200619-YBFB-e15-1'] = 0.00609
    J['200619-YBFB-e15-2'] = 0.00682
    J['200619-YBFB-e11-2'] = 0.00382
    J['200702-YBFB-Col-1'] = 1.62E-02  # - offset
    J['200703-YBFB-Col-2'] = 5.32E-02
    J['200703-YBFB-Col-3'] = 0.0151
    J['200703-YBFB-e15-3'] = 0.00634
    J['200703-YBFB-e11-4'] = 3.12E-02
    J['200703-YBFB-e11-3'] = 0.0035
    J['200703-YBFB-e11-2'] = 0.0028
    J['200721-YBFB-e15-5'] = 0.00403
    J['200721-YBFB-e11-2'] = 0.00203
    J['200721-YBFB-e15-2'] = 0.00438
    J['200721-YBFB-e15-1'] = 0.00266
    J['200724-YBFB-e11-2'] = 0.0042
    J['200724-YBFB-Col-1'] = 0.018
    J['200812-YBFB-e11-3'] = 1.10e-03
    J['200812-YBFB-e15-1'] = 0.00154
    J['200812-YBFB-e15-3'] = 0.00547

    _Jv['200219-YBFB-esb11'] = [0.00302, 0.00612, 0.00906]
    _Jv['160316#2'] = [3.74E-02, 3.94E-02, 4.42E-02, 4.66E-02, 4.93E-02, 5.08E-02, 6.30E-02, 9.01E-02]  # 300 kPa
    # _Jv['160316#2'] = [3.94E-02, 4.66E-02, 5.08E-02, 9.01E-02]  # un sur deux
    _Jv['2020-01-19-9h06'] = [0.005, 0.00567, 0.00601, 0.00633, 0.00664, 0.00698]
    # _Jv['2020-01-19-9h06'] = [0.00567, 0.00633, 0.00698] # un sur deux
    # _Jv['15012020-1045'] = [0.0065, 0.0073, 0.0091, 0.01, 0.0108, 0.0121, 0.01, 0.0155]
    _Jv['15012020-1045'] = [0.0065, 0.0073, 0.0091, 0.01, 0.0108, 0.0121, 0.0155]  # -1 pt
    _Jv['1045-sim'] = [0.00598743963242, 0.00615578009207, 0.00666497843767, 0.009325556553439999,
                                0.0107910981124, 0.0116187869283, 0.0133000529018]
    # _Jv['15012020-1045'] = [0.0065, 0.0073, 0.0091, 0.01, 0.0108, 0.0121] # -2 last pts
    _Jv['200219-YBFB-anac038'] = [0.0147, 0.0193, 0.0314]
    _Jv['200219-YBFB-Col'] = [0.00522, 0.00568, 0.00663, 0.00881]  # col + offset
    _Jv['200219-YBFB-horst1'] = [0.0121, 0.0176, 0.0234, 0.0335]
    _Jv['200619-YBFB-Col-1'] = [0.0167, 0.0168, 0.0189, 0.0198]
    # _Jv['200619-YBFB-Col-1'] = [0.0168, 0.0198]
    _Jv['200619-YBFB-Col-2'] = [0.0338, 0.0353, 0.038, 0.044]
    _Jv['200619-YBFB-e15-1'] = [0.00666, 0.0072, 0.0076, 0.00832]
    _Jv['200619-YBFB-e15-2'] = [0.00704, 0.00735, 0.00896, 0.00907]
    _Jv['200619-YBFB-e11-2'] = [0.00399, 0.00396, 0.00414]
    _Jv['200702-YBFB-Col-1'] = [1.74E-02, 1.93E-02, 2.11E-02, 2.72E-02]
    _Jv['200703-YBFB-Col-2'] = [0.0549, 0.0606, 0.0838]
    _Jv['200703-YBFB-Col-3'] = [0.017, 0.0199, 0.0222]
    _Jv['200703-YBFB-e15-3'] = [0.00642, 0.00681, 0.00735]
    _Jv['200703-YBFB-e11-4'] = [3.33E-02, 3.28E-02]
    _Jv['200703-YBFB-e11-3'] = [3.30E-03, 4.29E-03, 3.91E-03]
    _Jv['200703-YBFB-e11-2'] = [0.00257, 0.00249, 0.00274]
    _Jv['200721-YBFB-e15-5'] = [0.00442, 0.00437, 0.00436]
    _Jv['200721-YBFB-e11-2'] = [0.00226, 0.00246, 0.00242]
    _Jv['200721-YBFB-e15-2'] = [0.00451, 0.00537, 0.00997]
    _Jv['200721-YBFB-e15-1'] = [0.00294, 0.00312, 0.00369]
    _Jv['200724-YBFB-e11-2'] = [0.00471, 0.00478, 0.00518, 0.0053, 0.00564]
    _Jv['200724-YBFB-Col-1'] = [0.0205, 0.0243, 0.0277, 0.0347] # without the last pt
    _Jv['200812-YBFB-e11-3'] = [1.59e-03, 1.33e-03, 1.51e-03]
    _Jv['200812-YBFB-e15-1'] = [0.00193, 0.00203, 0.00256]
    _Jv['200812-YBFB-e15-3'] = [0.00558, 0.00575, 0.00609]

    cut_length['200219-YBFB-esb11'] = [0.0978, 0.0692, 0.0541]
    cut_length['160316#2'] = [0.13237, 0.12127, 0.09977, 0.08486, 0.07257, 0.0608, 0.052, 0.04252]
    # cut_length['160316#2'] = [0.12127, 0.08486, 0.0608, 0.04252]  # un sur deux
    cut_length['2020-01-19-9h06'] = [0.09124, 0.076704, 0.066162, 0.056812, 0.045564, 0.03707]
    # cut_length['2020-01-19-9h06'] = [0.076704, 0.056812, 0.03707]  # un sur deux
    # i['15012020-1045'] = [0.11059, 0.09159, 0.08232, 0.07078, 0.06082, 0.05338, 0.04282, 0.03649]
    cut_length['15012020-1045'] = [0.11059, 0.09159, 0.08232, 0.07078, 0.06082, 0.05338, 0.03649]  # -1 pt
    cut_length['1045-sim'] = [0.11059, 0.09159, 0.08232, 0.07078, 0.06082, 0.05338, 0.03649]
    # i['15012020-1045'] = [0.11059, 0.09159, 0.08232, 0.07078, 0.06082, 0.05338] # -2 last pts
    cut_length['200219-YBFB-anac038'] = [0.154, 0.117, 0.095]
    cut_length['200219-YBFB-Col'] = [0.1201, 0.0928, 0.0730, 0.0489]
    cut_length['200219-YBFB-horst1'] = [0.125, 0.106, 0.079, 0.046]
    cut_length['200619-YBFB-Col-1'] = [101.692e-3, 73.376e-3, 56.509e-3, 48.541e-3]
    # cut_length['200619-YBFB-Col-1'] = [73.376e-3, 48.541e-3]
    cut_length['200619-YBFB-Col-2'] = [112.3310e-3, 75.0720e-3, 60.5160e-3, 43.5200e-3]
    cut_length['200619-YBFB-e15-1'] = [116.809e-3, 85.251e-3, 62.335e-3, 42.325e-3]
    cut_length['200619-YBFB-e15-2'] = [96.991e-3, 70.135e-3, 48.791e-3, 36.769e-3]
    cut_length['200619-YBFB-e11-2'] = [0.104016, 0.068853, 0.04467]
    # i['200619-YBFB-e11-2'] = [0.068853, 0.04467] # - 1 pt
    cut_length['200702-YBFB-Col-1'] = [0.0978, 0.07496, 0.0547, 0.0333]
    cut_length['200703-YBFB-Col-2'] = [1.01E-01, 6.91E-02, 4.64E-02]
    cut_length['200703-YBFB-Col-3'] = [0.06533, 0.042995, 0.02946]
    cut_length['200703-YBFB-e15-3'] = [0.05426, 0.040365, 0.02971]
    cut_length['200703-YBFB-e11-4'] = [0.036322, 0.02888]
    cut_length['200703-YBFB-e11-3'] = [0.05117, 0.0327, 0.0222]
    cut_length['200703-YBFB-e11-2'] = [0.0546, 0.0408, 0.0298]
    cut_length['200721-YBFB-e15-5'] = [0.07068, 0.0462, 0.031]
    cut_length['200721-YBFB-e11-2'] = [0.080306, 0.040949, 0.02898]
    cut_length['200721-YBFB-e15-2'] = [0.06618, 0.05101, 0.039104]
    cut_length['200721-YBFB-e15-1'] = [.07099, 0.05413, 0.04169]
    cut_length['200724-YBFB-e11-2'] = [0.102385, 0.078216, 0.062506, 0.049686, 0.04045]
    cut_length['200724-YBFB-Col-1'] = [0.117942, 0.095206, 0.075436, 0.057301] # without the last pt
    cut_length['200812-YBFB-e11-3'] = [0.071707, 0.053612, 0.03859]
    cut_length['200812-YBFB-e15-1'] = [0.064336, 0.047067, 0.03706]
    cut_length['200812-YBFB-e15-3'] = [0.061457, 0.035634, 0.02968]

    if "esb11-random-" in filename:
        J[filename] = J['200219-YBFB-esb11']
        _Jv[filename] = _Jv['200219-YBFB-esb11']
        cut_length[filename] = cut_length['200219-YBFB-esb11']

    ### Maize
    J['4022020-841-maize'] = 0.157

    _Jv['4022020-841-maize'] = [0.1631, 0.1677, 0.1718, 0.1775, 0.2056, 0.35]

    cut_length['4022020-841-maize'] = [0.253, 0.203, 0.149, 0.118, 0.089, 0.04]

    for key in J.keys():
        if key in filename:
            plant = key
            break

    return J[plant], _Jv[plant], cut_length[plant]

def flux_for_k0_optim(filename = None, media = 'ctr-ctr'):
    J = {}
    P0 = {} # P @ Jv = 0 to account for osmotic effect on water potential
    media = media.lower()

    J['200219-YBFB-esb11'] = 0.00284
    J['200219-YBFB-esb11_old'] = 0.00284
    J['200219-YBFB-anac038'] = 0.0122
    J['200219-YBFB-horst1'] = 0.0116

    J['170426-full-archi-ch2D1'] = 1.96E-02
    J['170426-full-archi-ch2G1'] = 1.63E-02
    J['170426-full-archi-ch2G2'] = 1.41E-02
    J['170426-full-archi-ch3D1'] = 2.60E-02
    J['170426-full-archi-ch3G1'] = 1.97E-02
    J['170426-full-archi-ch3G2'] = 8.18E-03
    J['170426-full-archi-ch4D1b'] = 6.13E-03
    J['170426-full-archi-ch4G1c'] = 2.07E-02
    J['hydroroot-1ter'] = 4.52E-03
    J['160212-3-Col'] = 9.50E-03
    J['160212-4-Col'] = 9.88E-03
    J['hydroroot-5'] = 3.40E-03 # - offset 3.09E-03
    J['hydroroot-7'] = 1.38E-02 # - offset 1.55E-02
    J['hydroroot-8'] = 7.18E-03 # - offset 7.20E-03
    J['200812-YBFB-col-2'] = 2.96E-03
    J['200721-YBFB-Col-2'] = 1.18E-02

    J['180719-e15ch1-1'] = 4.21E-03
    J['180829-e15-3'] = 4.79E-03
    J['180829-e15-4'] = 5.27E-03

    J['150218-esk11-2'] = 1.42E-03
    J['150218-esk11-7'] = 8.64E-03
    J['220218-esk11-26'] = 1.47E-03
    J['220218-esk11-40'] = 1.71E-03
    J['220218-esk11-52'] = 8.55E-04
    J['220218-esk11-59'] = 1.84E-03
    J['180719-e11ch1-1'] = 2.14E-03
    J['180828-e11-1'] = 1.45E-03
    
    
    # Maize PR for primary root, sem for seminals
    # For semimals we may have Jv for 1 or 2 roots (left and right seminal)
    #       when we have 2 roots only the name  of the left one is given
    if media == 'ctr-ctr':
        # here I took the 1st measurement because some are done before ctr-150 measurements
        J['exp33-CTRpr1001-archi'] = 1.02E-01  # 1st measurement
        J['exp33-CTRpr3003-archi'] = 1.08E-01  # 1st measurement
        J['exp33-CTRpr5005-archi'] = 7.67E-02  # 1st measurement
        J['exp33-CTRpr8008-archi'] = 1.18E-01  # 1st measurement
        J['exp33-CTRpr10010-archi'] = 5.08E-02  # 1st measurement
        J['exp33-CTRpr7007-archi'] = 8.65E-02  # 1st measurement
        J['exp33-CTRpr2002-archi'] = 1.13E-01  # 1st measurement
        J['exp33-CTRpr4004-archi'] = 4.58E-02  # 1st measurement
        # J['exp33-CTRpr1001-archi'] = 1.06E-01  # measurement average
        # J['exp33-CTRpr3003-archi'] = 1.52E-01  # measurement average
        # J['exp33-CTRpr5005-archi'] = 6.86E-02  # measurement average
        # J['exp33-CTRpr8008-archi'] = 9.54E-02  # measurement average
        # J['exp33-CTRpr10010-archi'] = 4.57E-02  # measurement average
        # J['exp33-CTRpr1001-archi'] = 1.27E-01  # 1st measurement -P0
        # J['exp33-CTRpr3003-archi'] = 1.32E-01  # 1st measurement -P0
        # J['exp33-CTRpr5005-archi'] = 9.08E-02  # 1st measurement -P0
        # J['exp33-CTRpr8008-archi'] = 1.47E-01  # 1st measurement -P0
        # J['exp33-CTRpr10010-archi'] = 6.58E-02  # 1st measurement -P0
        # J['exp33-CTRpr7007-archi'] = 1.10E-01  # 1st measurement -P0
        
        J['exp33-CTRsem1001-L-archi'] = 3.59E-02  # 1st measurement
        J['exp33-CTRsem2002-L-archi'] = 1.52E-01  # 1st measurement
        J['exp33-CTRsem3003-L-archi'] = 2.62E-02  # 1st measurement
        J['exp33-CTRsem4004-L-archi'] = 3.75E-02  # 1st measurement
        J['exp33-CTRsem5005-L-archi'] = 9.54E-03  # 1st measurement
        J['exp33-CTRsem6006-L-archi'] = 2.80E-02  # 1st measurement
        J['exp33-CTRsem8008-L-archi'] = 2.72E-02  # 1st measurement
        J['exp33-CTRsem9009-L-archi'] = 1.51E-02 # 1st measurement

        # J['exp34-ctrpr1001-archi'] = 1.37E-01 # 1st measurement
        # J['exp34-ctrpr3003-archi'] = 2.01E-01 # 1st measurement
        # J['exp34-ctrpr6006-archi'] = 1.90E-01 # 1st measurement
        # J['exp34-ctrpr8008-archi'] = 1.07E-01 # 1st measurement
        # J['exp34-ctrpr5005-archi'] = 2.11E-01  # 1st measurement
        J['exp34-ctrpr1001-archi'] = 1.62E-01 # measurements average
        J['exp34-ctrpr3003-archi'] = 2.04E-01 # measurements average
        J['exp34-ctrpr6006-archi'] = 1.57E-01 # measurements average
        J['exp34-ctrpr8008-archi'] = 8.91E-02 # measurements average
        J['exp34-ctrpr5005-archi'] = 2.36E-01 # measurements average
        J['exp34-ctrpr2002-archi'] = 1.22E-01 # 1st measurement
        J['exp34-ctrpr9009-archi'] = 3.17E-01 # 1st measurement
        J['exp34-ctrpr47004-archi'] = 2.47E-01 # 1st measurement
        
        J['exp34-ctrsem1001-L-archi'] = 2.03E-02 # 1st measurement
        J['exp34-ctrsem2002-L-archi'] = 4.46E-02 # 1st measurement
        J['exp34-ctrsem3003-archi'] = 3.10E-02 # 1st measurement
        J['exp34-ctrsem6005-L-archi'] = 7.26E-02 # 1st measurement
        J['exp34-ctrsem7006-archi'] = 5.02E-02 # 1st measurement
        J['exp34-ctrsem8007-L-archi'] = 2.28E-01 # 1st measurement
        J['exp34-ctrsem9008-L-archi'] = 9.54E-02 # 1st measurement
        J['exp34-ctrsem10009-archi'] = 5.99E-02 # 1st measurement
        J['exp34-ctrsem11010-archi'] = 2.73E-02 # 1st measurement
        J['exp34-ctrsem12011-l-archi'] = 1.06E-01 # 1st measurement

    elif media == 'ctr-150':
        # here I took the last measurement because of the dynamical aspect
        J['exp33-CTRpr2002-archi'] = 3.73E-02 # last point only
        J['exp33-CTRpr4004-archi'] = 1.70E-02 # last point only
        J['exp33-CTRpr6006-archi'] = 2.95E-02 # last point only
        J['exp33-CTRpr9009-archi'] = 1.42E-02 # last point only
        P0['exp33-CTRpr2002-archi'] = 77.21 # last point only
        P0['exp33-CTRpr4004-archi'] = -10.05 # last point only
        P0['exp33-CTRpr6006-archi'] = 19.02 # last point only
        P0['exp33-CTRpr9009-archi'] = -63.44 # last point only
        # J['exp33-CTRpr2002-archi'] = 3.32E-02  # average whithout 1st point
        # J['exp33-CTRpr4004-archi'] = 1.75E-02 # average without 1st point
        # J['exp33-CTRpr6006-archi'] = 2.73E-02 # average without 1st point
        # J['exp33-CTRpr9009-archi'] = 1.38E-02 # average without 1st point
        # P0['exp33-CTRpr2002-archi'] = 39.38 # average without 1st point
        # P0['exp33-CTRpr4004-archi'] = -14.3 # average without 1st point
        # P0['exp33-CTRpr6006-archi'] = -25.48 # average without 1st point
        # P0['exp33-CTRpr9009-archi'] = -56.45 # average without 1st point

        J['exp33-CTRsem2002-L-archi'] = 2.27E+00
        J['exp33-CTRsem4004-L-archi'] = 1.19E-02
        J['exp33-CTRsem6006-L-archi'] = 1.89E-02
        J['exp33-CTRsem8008-L-archi'] = 1.52E-02
        P0['exp33-CTRsem2002-L-archi'] = 62.5
        P0['exp33-CTRsem4004-L-archi'] = -1.81
        P0['exp33-CTRsem6006-L-archi'] = -132.96
        P0['exp33-CTRsem8008-L-archi'] = -180.17

        J['exp34-ctrpr2002-archi'] = 3.22E-02
        J['exp34-ctrpr47004-archi'] = 2.96E-02
        J['exp34-ctrpr9009-archi'] = 7.28E-02
        P0['exp34-ctrpr2002-archi'] = -17.61
        P0['exp34-ctrpr47004-archi'] = 28.41
        P0['exp34-ctrpr9009-archi'] = -103.46
        
        J['exp34-ctrsem2002-L-archi'] = 1.27E-02
        J['exp34-ctrsem6005-R-archi'] = 2.25E-02
        J['exp34-ctrsem8007-L-archi'] = 1.56E-02
        P0['exp34-ctrsem2002-L-archi'] = -108.11
        P0['exp34-ctrsem6005-R-archi'] = 23.18
        P0['exp34-ctrsem8007-L-archi'] = -86.45
        P0['exp34-ctrsem10009-archi'] = -30.9 # Last measurement
        P0['exp34-ctrsem11010-archi'] = -1092.72 # Last measurement
        P0['exp34-ctrsem12011-l-archi'] = -221.76 # Last measurement
        J['exp34-ctrsem10009-archi'] = 8.47E-03 # Last measurement
        J['exp34-ctrsem11010-archi'] = 5.66E-03 # Last measurement
        J['exp34-ctrsem12011-l-archi'] = 2.43E-02 # Last measurement

    elif media == '150-150':
        # here I took the 1st measurement because some are done before 150-ctr measurements
        J['exp33-150PR1001-archi'] = 1.08E-02
        J['exp33-150PR2003-archi'] = 6.67E-03
        J['exp33-150PR3004-archi'] = 2.05E-02
        J['exp33-150PR4005-archi'] = 1.14E-02
        J['exp33-150PR5006-archi'] = 1.13E-02
        J['exp33-150PR6009-archi'] = 4.82E-03
        P0['exp33-150PR1001-archi'] = 22.03
        P0['exp33-150PR2003-archi'] = 261.97
        P0['exp33-150PR3004-archi'] = 77.41
        P0['exp33-150PR4005-archi'] = -105.9
        P0['exp33-150PR5006-archi'] = -3.79
        P0['exp33-150PR6009-archi'] = -48.87

        J['exp33-150sem3003-archi'] = 1.02E-02
        J['exp33-150sem4004-L-archi'] = 5.17E-03
        J['exp33-150sem5005-L-archi'] = 4.25E-03
        J['exp33-150sem6006-L-archi'] = 1.34E-02
        P0['exp33-150sem3003-archi'] = -70.9
        P0['exp33-150sem4004-L-archi'] = -225.76
        P0['exp33-150sem5005-L-archi'] = -475.7
        P0['exp33-150sem6006-L-archi'] = -177.95

        J['exp34-150pr1'] = 1.30E-02 # 1st measurement
        J['exp34-150pr2'] = 2.38E-02 # 1st measurement
        J['exp34-150pr3'] = 2.37E-02 # 1st measurement
        J['exp34-150pr5'] = 1.88E-02 # 1st measurement
        J['exp34-150pr6'] = 5.79E-02 # 1st measurement
        J['exp34-150pr7'] = 1.56E-02 # 1st measurement
        J['exp34-150pr8'] = 1.65E-02 # 1st measurement
        J['exp34-150pr9'] = 2.45E-02 # 1st measurement
        J['exp34-150pr10'] = 1.85E-02 # 1st measurement
        J['exp34-150pr11'] = 1.46E-02 # 1st measurement
        J['exp34-150pr13'] = 1.23E-02 # 1st measurement
        J['exp34-150pr14'] = 3.13E-02 # 1st measurement
        J['exp34-150pr15'] = 5.17E-02 # 1st measurement

        P0['exp34-150pr1'] = -9.08E+01 # 1st measurement
        P0['exp34-150pr3'] = -2.93E+01 # 1st measurement
        P0['exp34-150pr5'] = -1.13E+02 # 1st measurement
        P0['exp34-150pr7'] = -1.22E+02 # 1st measurement
        P0['exp34-150pr9'] = -1.05E+02 # 1st measurement
        P0['exp34-150pr13'] = -3.77E+02 # 1st measurement
        P0['exp34-150pr10'] = 7.88E+01 # 1st measurement
        P0['exp34-150pr11'] = -1.42E+02 # 1st measurement
        P0['exp34-150pr2'] = 1.08E+02 # 1st measurement
        P0['exp34-150pr6'] = -4.79E+01 # 1st measurement
        P0['exp34-150pr8'] = 9.59E+01 # 1st measurement
        P0['exp34-150pr14'] = -5.37E+01 # 1st measurement
        P0['exp34-150pr15'] = -3.66E+00 # 1st measurement

        
        J['exp34-150sem3003-archi'] = 3.12E-03
        J['exp34-150sem6005-archi'] = 3.59E-03
        J['exp34-150sem7006-L-archi'] = 1.80E-02
        J['exp34-150sem8007-archi'] = 6.73E-03
        J['exp34-150sem9008-archi'] = 2.45E-02
        J['exp34-150sem11010-archi-L'] = 6.85E-02
        P0['exp34-150sem3003-archi'] = -252.42
        P0['exp34-150sem6005-archi'] = -90.77
        P0['exp34-150sem7006-L-archi'] = -8.63
        P0['exp34-150sem8007-archi'] = -112.53
        P0['exp34-150sem9008-archi'] = -36.68
        P0['exp34-150sem11010-archi-L'] = 67.41


    elif media == '150-ctr':
        # here I took the last measurement because of the dynamical aspect
        J['exp33-150PR2003-archi'] = 1.86E-02
        J['exp33-150PR4005-archi'] = 3.58E-02
        J['exp33-150PR6009-archi'] = 2.14E-02
        P0['exp33-150PR2003-archi'] = -177.25
        P0['exp33-150PR4005-archi'] = -37.6
        P0['exp33-150PR6009-archi'] = -3.15

        J['exp33-150sem2002-archi'] = 2.57E-03
        J['exp33-150sem4004-L-archi'] = 1.15E-02
        J['exp33-150sem6006-L-archi'] = 1.18E-02
        P0['exp33-150sem2002-archi'] = -218.76
        P0['exp33-150sem4004-L-archi'] = -204.36
        P0['exp33-150sem6006-L-archi'] = -153.73

        J['exp34-150pr6005-archi'] = 6.82E-02
        J['exp34-150pr8007-archi'] = 2.80E-02
        J['exp34-150pr11010-archi'] = 8.24E-02
        P0['exp34-150pr6005-archi'] = -3.66
        P0['exp34-150pr8007-archi'] = -342.09
        P0['exp34-150pr11010-archi'] = 23.94

        J['exp34-150sem6005-archi'] = 1.36E-02
        J['exp34-150sem8007-archi'] = 6.57E-03
        J['exp34-150sem11010-archi-L'] = 6.53E-03
        P0['exp34-150sem6005-archi'] = 22.98
        P0['exp34-150sem8007-archi'] = -122.3
        P0['exp34-150sem11010-archi-L'] = -145.41

    elif media == 'ctr-ctr-pr-266-835':
        # simulations done with k0=266 for PR and k0=835.24 for LR
        # Kx [0, 0.030196, 0.046758, 0.067525, 0.3] # remark : last point at 0.3 !!!
        #    [0.00132, 0.00129,0.00127,0.04230,0.04207]
        J['exp34-ctrpr1001-archi'] = 0.15339471092
        J['exp34-ctrpr2002-archi'] = 0.170172028187
        J['exp34-ctrpr3003-archi'] = 0.179452306987
        J['exp34-ctrpr6006-archi'] = 0.201299886492
        J['exp34-ctrpr8008-archi'] = 0.159410015208
        J['exp34-ctrpr9009-archi'] = 0.212967699241
        J['exp34-ctrpr47004-archi'] = 0.157448316912
    
    elif media == 'sim-455-87.5':
        # simulations done with k0=87.5 for PR and k0=455 for LR
        # Kx [0, 0.030196, 0.046758, 0.067525, 0.5]
        #    [0.00132, 0.00129,0.00127,0.04230,0.04207]
        # generated done with the experimental primary root length
        J['generated-Exp33-16108107'] = 0.094
        J['generated-Exp33-28601240'] = 0.103
        J['generated-Exp33-40025372'] = 0.092
        J['generated-Exp33-44487771'] = 0.109
        J['generated-Exp33-65754910'] = 0.086
        J['generated-Exp33-77288684'] = 0.091
        J['generated-Exp33-86717948'] = 0.102
        J['generated-Exp33-96647360'] = 0.111
        J['generated-Exp33-96768222'] = 0.083
        J['generated-Exp34-95661757.csv'] = 0.1344
        J['generated-Exp34-76244324.csv'] = 0.1289
        J['generated-Exp34-48438801.csv'] = 0.1061
        J['generated-Exp34-79470944.csv'] = 0.1672
        J['generated-Exp34-96810412.csv'] = 0.1539
        J['generated-Exp34-72001363.csv'] = 0.1593
        J['generated-Exp34-78987615.csv'] = 0.1584

        J['exp34-ctrpr1001-archi'] = 1.19E-01
        J['exp34-ctrpr2002-archi'] = 1.19E-01
        J['exp34-ctrpr3003-archi'] = 1.34E-01
        J['exp34-ctrpr6006-archi'] = 1.55E-01
        J['exp34-ctrpr8008-archi'] = 1.22E-01
        J['exp34-ctrpr9009-archi'] = 1.66E-01
        J['exp34-ctrpr47004-archi'] = 1.20E-01

    Jv = None
    PJv0 = 0.0
    for key in J.keys():
        if key in filename:
            plant = key
            Jv = J[key]
            if key in P0.keys(): PJv0 = P0[key] * 1e-3
            break

    return Jv, PJv0

