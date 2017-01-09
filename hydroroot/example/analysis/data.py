""" Define data from measurement and experiment.

Provides methods for the different genotypes.
"""

def col():
    l={}

    l['median'] = ([x/100. for x in range(0,101,5)],
                [0.0, 0., 0., 0.00032, 0.000835, 0.00076, 0.002665, 0.003825, 0.006775,
                 0.00907, 0.01282, 0.01193, 0.01129, 0.010465, 0.015035, 0.01016, 0.016275,
                 0.01536, 0.05308, 0.07251,0.07971])
    l['low']= ([x/100. for x in range(0,101,5)],
                [0.0, 0., 0., 0.00032, 0.00032, 0.000515, 0.002252, 0.001312, 0.003844,
                 0.006908, 0.005966, 0.00089, 0.001751, 0.001921, 0.002584, 0.003325, 0.005306,
                 0.009524, 0.046133, 0.056527,0.074152])
    l['high']= ([x/100. for x in range(0,101,5)],
                [0.0, 0., 0., 0.00032, 0.001465, 0.0012, 0.003725, 0.007105, 0.01336,
                 0.009486, 0.02051, 0.022757, 0.023004, 0.027506, 0.0215, 0.02938, 0.039895,
                 0.050636, 0.062388, 0.07646,0.09448])
    return l

def primary_length_law(length, genotype='col', type='median'):
    #print 'Length law at ', length
    #length=0.12
    if genotype == 'col':
        l = col()
    else:
        raise Exception('Length law is not defined for the genotype %s'%genotype)

    x, y = l[type]

    return [xx*length for xx in x], y

def higher_order_length_law(max_length, genotype='col', type='median'):
    return primary_length_law(max_length, genotype, type)

# radial
def radial(v=300):
    xr = [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16]
    yr = [v]*len(xr)
    return xr, yr

def axial_law(genotype='col'):
    # Remove last value
    a={}
    a['col'] = (
        [0., 0.03,  0.06, 0.09, 0.12, 0.15, 0.18],
        [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4,292.6e-4,262.5e-4,511.1e-4]
    )

    a['esk11'] = (
        [0., 0.03,  0.06, 0.09, 0.12],
        [2.3e-4, 15.9e-4, 158.9e-4, 216.4e-4,192.e-4]
    )

    a['esk15'] = (
        [0., 0.03,  0.06, 0.09, 0.12, 0.15],
        [0.7e-4, 3.9e-4, 168.8e-4, 220.5e-4,239.3e-4, 219.5e-4]
    )

    a['irx34'] = (
        [0., 0.03,  0.06, 0.09, 0.12],
        [5.2e-4, 11.4e-4, 19.8e-4, 119.6e-4,359.9e-4]
    )

    a['pip2122'] = a['col']

    return a[genotype]
"""
r = {}
r['col'] = radial(300)
r['pip2122'] = radial(239)
r['esk11'] = radial(373)
r['esk15'] = radial(518)
r['irx34'] = radial(300)

# axial

"""
