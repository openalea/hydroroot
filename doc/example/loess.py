import copy
import glob

import pandas as pd
import numpy as np
from scipy import stats
import scipy.interpolate
import statsmodels.api as sm
import matplotlib.pyplot as plt
import math
# import agramfort_lowess as al
# import lowess as duhs
# import skmisc, math
# from skmisc.loess import loess

"""
The confidence interval:
-	Ref WEISBERG, “Applied Linear Regression” 3d ed., Wiley, pp 277-278 1997 to evaluate the std dev
	Smoothing done on ri ² = (yi – yfiti)² giving the smoothed estimate of the variance, which gives the estimate of the standard deviation si
-	Then the interval is given by ±t*si/sqrt(n) t from the 5% cut on the t distribution
"""

print('***************************')
print('*** check if at the base some LR are > than lcut')
print('***************************')

lowess = sm.nonparametric.lowess
directory_name='/home/bauget/Documents/Rice/Experiments/Gaelle/results/'
# directory_name=''

# output_name= 'K-241220-CL.csv'
# output_name = directory_name + 'K-lowess-optim-J_P_nonuniform-150.csv'
# output_name = 'lowess-PEG-CTR-Jv.csv'

# x_col_name = "x pr"
# y_col_name = "K pr"
# x_col_name = "x" # no-solute script
# y_col_name = "K optimized" # no-solute script
x_col_name = "x" # no-solute script
y_col_name = "K optimized" # no-solute script
# x_col_name = "dp"
# y_col_name = "Jexp Lpr"
# kr_col_name = 'kpr'
kr_col_name = 'k (10-9 m/s/MPa)'
_mean = {}
for key in ['exp', 'kr', 'k SE', 'Js', 'Js SE', 'Ps', 'Ps SE']:
    _mean[key] = []

frac = 2.0/3.0
iteration = 2

exclude_list = ['results-exp03-k2.csv']
_list = ['*k*.csv']
for j in _list :
    print(j)
    kr = []
    js = []
    ps = []
    files_list=glob.glob(directory_name + j)
    # output_name = j.replace('*','-').replace('tomato', 'K')
    output_name = j.replace('*','-').replace('solute', 'K')
    dk2 = None
    # output_name = directory_name + 'K-lowess-line-' + a + '-' + str(j) + '.csv'
    for i in range(len(files_list)):
        f_name=files_list[i]
        if f_name.replace(directory_name, '') in exclude_list: continue
        dk1 = pd.read_csv(f_name, sep = ',')
        kr.append(dk1[kr_col_name][0])
        # js.append(dk1['Js'][0])
        # ps.append(dk1['Ps'][0])
        if dk1[x_col_name].iloc[-1] == 0.0:
            n = len(dk1[x_col_name])
            dk1=dk1.drop([n-1])
        # normalization
        # xmin = dk1[x_col_name].min()
        # xmax = dk1[x_col_name].max()
        # ymin = dk1[y_col_name].min()
        # ymax = dk1[y_col_name].max()
        #
        # dnorm=(dk1[y_col_name]-ymin)/(ymax-ymin)
        # dk1[y_col_name] = dnorm
        # dnorm=(dk1[x_col_name]-xmin)/(xmax-xmin)
        # dk1[x_col_name] = dnorm

        dk2 = pd.concat([dk1,dk2])


    dksorted = dk2[[x_col_name,y_col_name]]
    dksorted = dksorted.dropna(subset = [x_col_name]) #, inplace=True)
    dksorted = dksorted.sort_values(x_col_name)

    x = np.array(dksorted[x_col_name])
    y = np.array(dksorted[y_col_name])
    n = np.size(y)



    # for small samples (<50) we use t-statistics
    # for 95% confidence interval, alpha = 5% = 0.05 and alpha/2 = 0.025
    alpha = 0.05
    t = stats.t.ppf(1 - alpha/2, n - 2)
    fig = plt.figure()
    plt.title(output_name)

    z1 = lowess(y,x,frac=frac,it=iteration,is_sorted=True,return_sorted=True) #, delta = 0.01)
    res = (y - z1[:,1])**2.0
    res_est = lowess(res, x, frac=frac, it=iteration, is_sorted=True, return_sorted=True)
    std = np.sqrt(np.abs(t*res_est/math.sqrt(n))) # +/- conf interval
    std[:,0] = x

    # derivative
    z2 = np.unique(z1, axis = 0)
    der = np.diff(z2[:,1])/np.diff(z2[:,0])

    xlo = np.unique(z1, axis = 0)[:, 0]
    ylo = np.unique(z1, axis = 0)[:, 1]
    intpl = scipy.interpolate.UnivariateSpline(xlo,ylo)
    # rbf = scipy.interpolate.Rbf(np.unique(z1, axis = 0)[:, 0], np.unique(z1, axis = 0)[:, 1])
    dp=0.025
    xlo=np.array(np.arange(0.0, 1.0, dp))
    ylo = intpl(xlo)
    dx= np.gradient(xlo)
    dydx = np.gradient(ylo) / dx
    d2ydx = np.gradient(dydx) / dx
    curvature = d2ydx / ((1.0 + dydx ** 2.0) ** (3.0 / 2.0))
    R = 1/curvature

    plt.plot(x, y, 'o')
    plt.plot(z1[:,0],z1[:,1])
    # plt.plot(x,z1)
    plt.fill_between(z1[:,0], z1[:,1]-std[:,1], z1[:,1]+std[:,1], alpha = .33)
    # plt.fill_between(eval_x, bottom, top, alpha = 0.5, color = "b")
    fig.patch.set_facecolor('lightgrey')
    fig.tight_layout()
    # fig2 = plt.figure()
    # plt.plot(xlo, curvature)
    # fig2.patch.set_facecolor('lightgrey')
    # fig2.tight_layout()

    # dout = pd.DataFrame(np.hstack([z1,std]))
    # dout = dout.drop_duplicates()
    dout = pd.DataFrame()
    dout['x'] = z1[:,0]
    dout['K'] = z1[:,1]
    dout['-std'] = z1[:,1]-std[:,1]
    dout['+std'] = z1[:,1]+std[:,1]
    dout = dout.drop_duplicates()
    if output_name is not None: dout.to_csv(output_name, index=False)
    # dkr = pd.DataFrame(kr)
    # djs = pd.DataFrame(js)
    # dps = pd.DataFrame(ps)
    # d2 = pd.concat([dkr, djs, dps], axis = 1)
    # d2.to_csv('out2.csv')
    array = np.array(kr)
    _mean['kr'].append(array.mean())
    _mean['k SE'].append(array.std(ddof=1)/np.sqrt(len(array)))
    array = np.array(js)
    _mean['Js'].append(array.mean())
    _mean['Js SE'].append(array.std(ddof=1)/np.sqrt(len(array)))
    array = np.array(ps)
    _mean['Ps'].append(array.mean())
    _mean['Ps SE'].append(array.std(ddof=1)/np.sqrt(len(array)))
    _mean['exp'].append(output_name.replace(directory_name, '').replace('.csv', '(n=' + str(len(array)) + ')'))

d=pd.DataFrame(_mean, columns = ['exp', 'kr', 'k SE', 'Js', 'Js SE', 'Ps', 'Ps SE'])
d.to_csv('out.csv', index = False)
