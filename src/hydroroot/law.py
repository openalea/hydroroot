# -*- coding: utf-8 -*-
"""
Define a length law of LRs using equal amplitudes method from measured data.

Created on Tue Feb 17 12:56:51 2015

@author: ndour
"""

import random

import numpy as np
import pylab

from hydroroot import length


def expovariate_law(data_xy, size=5e-2, scale_x=1e-2, scale_y=1e3, plot=False):
    """
    Fit a spline law from measured data by adding stochasticity.

    To compute the law, data are first sampled using equal amplitude method. Then, a mean is computed on each sample.
    Then, an expovariate distribution is simulated from the mean of each sample.
    Finally a spline is interpolated based on the simulated data.

    Parameters
    ==========
        - data_xy:
        - size:
        - scale_x:
        - scale_y:
        - plot:

    Example
    =======

        >>> filename = 'lr_length_law_data.csv'
        >>> xy = readCSVFile(filename)
        >>> law = fit_law(data_xy=xy, size=5e-2)
    """

    # sort by position
    xy = data_xy
    xy.sort(axis=0)
    xy = xy.tolist()
    x, y = zip(*xy)  # Separate x and y coordiantes
    x = np.array(x) * scale_x
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    X, values = discretize(x, y, size=size)

    Y = [np.mean(ys) for ys in values]
    YY = [(random.expovariate(1. / v) if v > 0 else 0.) for v in Y]

    if plot:
        Y_max = [max(ys) for ys in values]
        Y_min = [min(ys) for ys in values]

        pylab.plot(x, y)
        pylab.plot(X, Y_max)
        pylab.plot(X, Y_min)
        pylab.plot(X, YY)

    law = length.fit_law(X, YY)

    return law


def discretize(x, y, size=5e-2):
    """ Discretize by 5cm-intervals by using equal amplitudes method"""

    m, M = min(x), max(x)

    delta = (M - m) / float(size)
    points = np.linspace(m, M, delta).tolist()

    #points = [(m + i * size) for i in range(1, nb_class - 1)]

    #points.insert(0, m)
    #points.append(M)
    nb_points = len(points)
    intervals = [(points[i], points[i + 1]) for i in range(nb_points-1)]

    ys = [[y[i] for i, p in enumerate(x) if p1 <= p <= p2] for p1,p2 in intervals]

    ys.insert(0, [0.])


    zz = [(points[i], y) for i, y in enumerate(ys) if y]

    return zip(*zz)


def multi_law(x, y, size=5e-2, scale_x=0.16/100., scale_y=1e-3, plot=False):
    """

    """
    x = np.array(x) * scale_x
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    X, values = discretize(x, y, size)
    Y = [np.mean(ys) for ys in values]
    YY = [(random.expovariate(1. / v) if v > 0 else 0.) for v in Y]


    Y_max = [max(ys) for ys in values]
    Y_min = [min(ys) for ys in values]

    if plot:
        #pylab.plot(x, y, label='data')
        pylab.plot(X, Y_max, label='max')
        pylab.plot(X, Y_min, label='min')
        pylab.plot(X, YY, label='mean')
        pylab.legend()

    return (X, Y_min), (X, Y_max), (X,YY)


def histo_relative_law(x, y, size=5e-2, scale_x=1., scale_y=1e-3, scale=1e-4, plot=False, uniform=False):
    """ Return a length law from [0,1] to absolute length.

    Algorithm:
      - First, discretize the X values in different intervals of size `size`.
      - Compute the histogram from the set of points include in each interval.
      - Return a function that compute a value in a given histogram
    """

    x = np.array(x) * scale_x
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    X, values = discretize(x, y, size)

    means = [np.mean(ys) for ys in values]

    def return_law(position, scale=scale):
        for i, x_min in enumerate(X):
            if position*scale <= x_min:
                break
        index = max(i-1, 0)
        points = values[index]
        n = len(points)

        if n == 0:
            return 0.

        if uniform=='expo':
            v = means[index]
            length = random.expovariate(1. / v) if v > 0 else 0.
            # shoud not exceed the law, some randomness around length is done in markov, branching_variability
            if length > max(points):
                length = max(points)
        elif not uniform:
            index_value = random.randint(0,n-1)
            length = points[index_value]
        else:
            min_y = min(points)
            max_y = max(points)
            length = min_y + (max_y-min_y)*random.random()

        return length/scale

    return return_law

def reference_relative_law(x, y, size=5e-2, scale_x=1., scale_y=1e-3):
    """ Return a length law from [0,1] to absolute length.

    Algorithm:
      - First, discretize the X values in different intervals of size `size`.
      - Compute the histogram from the set of points include in each interval.
      - Return a function that compute a value in a given histogram
    """

    x = np.array(x) * scale_x
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    X, values = discretize(x, y, size)

    means = [np.mean(ys) for ys in values]

    return length.fit_law(X, means, ext=2)


