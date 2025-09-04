# -*- coding: utf-8 -*-
"""
Define a length law of LRs using equal amplitudes method from measured data.

Created on Tue Feb 17 12:56:51 2015

@author: ndour
"""

import random

import numpy as np
import pylab

from openalea.hydroroot import length


def expovariate_law(data_xy, size=5e-2, scale_x=1e-2, scale_y=1e3, plot=False):
    """Fit a spline law from measured data by adding stochasticity.
    
    To compute the law, data are first sampled using equal amplitude method. Then, a mean is computed on each sample.
    Then, an expovariate distribution is simulated from the mean of each sample.
    Finally a spline is interpolated based on the simulated data.

    :param data_xy: the data to fit
    :param size: the sample size (Default value = 5e-2)
    :param scale_x: a scaling factor on x array (Default value = 1e-2)
    :param scale_y: a scaling factor on x array (Default value = 1e3)
    :param plot: True plot the law (Default value = False)
    :param Example: 

    >>> filename = 'lr_length_law_data.csv'
        >>> xy = readCSVFile(filename)
        >>> law = fit_law(data_xy=xy, size=5e-2)
    """

    # sort by position
    xy = data_xy
    xy.sort(axis=0)
    xy = xy.tolist()
    x, y = list(zip(*xy))  # Separate x and y coordiantes
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
    """Discretize by intervals of size `size` by using equal amplitudes method

    :param x: (float list) - abscissa
    :param y: (float list)
    :param size: bins size (Default value = 5e-2)

    """

    m, M = min(x), max(x)

    # F. Bauget 2022-03-15: python 2 to 3
    #   - before in numpy.linspace the 3ed arguments was tranformed to int in the numpy routine
    #   - now should be done before call
    # delta = (M - m) / float(size)
    delta = int((M - m) / float(size))
    points = np.linspace(m, M, delta).tolist()

    #points = [(m + i * size) for i in range(1, nb_class - 1)]

    #points.insert(0, m)
    #points.append(M)
    nb_points = len(points)
    intervals = [(points[i], points[i + 1]) for i in range(nb_points-1)]

    ys = [[y[i] for i, p in enumerate(x) if p1 <= p <= p2] for p1,p2 in intervals]

    ys.insert(0, [0.])


    zz = [(points[i], y) for i, y in enumerate(ys) if y]

    return list(zip(*zz))


def multi_law(x, y, size=5e-2, scale_x=0.16/100., scale_y=1e-3, plot=False):
    """
    deprecated
    :param x: 
    :param y: 
    :param size: the bin size (Default value = 5e-2)
    :param scale_x:  (Default value = 0.16/100.)
    :param scale_y:  (Default value = 1e-3)
    :param plot:  (Default value = False)

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
    """
    Return a function return_law(position, scale=scale) that compute a value from y values (see below) binned according to x divided in bins of size `size`.

    :param x: (list of float)
    :param y: (list of float)
    :param size: (float) - bins size (Default value = 5e-2)
    :param scale_x: (float) - a factor that multiplies x values (Default value = 1.)
    :param scale_y: (float) - a factor that multiplies y values (Default value = 1e-3)
    :param scale: (float) - a number that divides the length given by the function (Default value = 1e-4)
    :param plot: unused (Default value = False)
    :param uniform: (string or boolean) - 'expo', True or False see below (Default value = False)
    :returns:
        - return_law

    Algorithm:
      - First, discretize the x values in different intervals of size `size`.
      - Compute the mean of the set of points included in each interval.
      - Return a function that computes for a given interval at position `position` :
        - if uniform='expo', randomly a value from an exponential distribution with mean equals to <y> in this interval
        - if uniform=False, randomly one of the y values in the interval
        - if uniform=True, a value randomly chosen between min(y) and max(y) in the interval

    """

    x = np.array(x) * scale_x
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    X, values = discretize(x, y, size)

    means = [np.mean(ys) for ys in values]

    def return_law(position, scale=scale):
        """

        :param position: (float)
        :param scale: - (Default value = scale)

        """
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
    """
    Return a spline that interpolates the binned mean of y, see below for the calculation.

    :param x: (list of float)
    :param y: (list of float)
    :param size: (float) - bins size (Default value = 5e-2)
    :param scale_x: (float) - a factor that multiplies x values (Default value = 1.)
    :param scale_y: (float) - a factor that multiplies y values (Default value = 1e-3)
    :returns:
        - spline

    :Algorithm:
      - First, discretize the X values in different intervals of size `size`.
      - Compute the mean of the set of points included in each interval.
      - interpolate with a spline y(x)

    """

    x = np.array(x) * scale_x
    y = np.array(y) * scale_y
    x = x.tolist()
    y = y.tolist()

    X, values = discretize(x, y, size)

    means = [np.mean(ys) for ys in values]

    return length.fit_law(X, means, ext=2)


def length_law(pd, scale_x = 1 / 100., scale_y = 1., scale = 1e-4, uniform = 'expo', size = 5):
    """Build the function giving the lateral length according to its position on the parent branch

    :param pd: DataFrame
    :param scale_x: (float) - a factor that multiplies size and x values (Default value = 1 / 100.)
    :param scale_y:  (float) - a factor that multiplies y values (Default value = 1.)
    :param scale: (float) number that will divide the given length, so for instance it used to
     divide the returned length (by `histo_relative_law`) by the segment length to get it in number of segment (Default value = 1e-4)
    :param uniform: boolean or string (Default value = 'expo') see :func:`~law.py.histo_relative_law`
    :param size:  (Default value = 5)
    :returns:
        - a function giving the lateral length according to its position

    Remark:
    This is specific to the length law files:
         * 1st col: "LR_length_mm(mm)" lateral lengths in mm
         * 2nd col: "relative_distance_to_tip" relative distance to tip in % so between 0 and 100.

    """
    x = pd.relative_distance_to_tip.tolist()
    y = pd.LR_length_mm.tolist()

    # size of the windows: in %
    size *= scale_x
    # TODO : change '1.e-3 * scale_y' to scale_y
    _length_law = histo_relative_law(x, y,
                                     size = size,
                                     scale_x = scale_x,
                                     scale_y = 1.e-3 * scale_y,
                                     scale = scale,
                                     plot = False,
                                     uniform = uniform)
    return _length_law
