__author__ = "Jwely"

import numpy as np
import math


def dbz(a, b):
    """
    divides two nd arrays element wise, but returns 0 if either denomenator
    or numerator is zero.
    :param a:  the numerator
    :param b:  the denominator
    :return:   a / b
    """

    if isinstance(a, np.ma.MaskedArray):
        mask = a.mask
        a = a.data
    if isinstance(b, np.ma.MaskedArray):
        b = b.data

    q = a / b
    q[q > 1e100] = 0
    q[q == np.nan] = 0
    q[q == -np.inf] = 0
    q[q == np.inf] = 0
    return np.ma.masked_array(q, mask=mask)


def get_spatial_derivative(f, x, y):
    """
    Returns a derivative of f with respect to x and y for discreet nd arrays of points.
    :param f:   2d array of some function f values
    :param x:   2d array of x coordinate values
    :param y:   2d array of y coordinate values
    """

    mask = np.zeros(f.shape)
    if isinstance(f, np.ma.MaskedArray):
        mask = f.mask

    # derivatives in terms of index locations i,j
    dfdi, dfdj = np.gradient(f)
    dxdi, dxdj = np.gradient(x)
    dydi, dydj = np.gradient(y)

    dfdx = dbz(dfdj, dxdj) + dbz(dfdi, dxdi)
    dfdy = dbz(dfdi, dydi) + dbz(dfdj, dydj)

    return dfdx, dfdy


# testing area
if __name__ == "__main__":


    # test out both cylindrical and cartesian
    xl = np.linspace(1, 100, 100)
    yl = np.linspace(1, 100, 100)
    x, y = np.meshgrid(xl, yl)
    xi, yi = (1.1, 1.1)
    x -= xi
    y -= yi

    r = ((x) ** 2 + (y) ** 2) ** 0.5
    t = np.arctan2(y, x)

    f = 5 * r + 6 * t
    dfdx, dfdy = get_spatial_derivative(f, x, y)
    dfdr = (dbz(dfdy, np.sin(t)) + dbz(dfdx, np.cos(t))) / 2
    dfdt = (dbz(dfdy, (1/r) * np.cos(t)) + dbz(dfdx, (-1/r) * np.sin(t))) / 2
    print np.percentile(dfdr, 50)
    print np.percentile(dfdt, 50)

