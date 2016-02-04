__author__ = "Jwely"

import numpy as np


def divide_by_zero(a, b, mask):
    """
    divides two nd arrays element wise, but returns 0 if either denomenator
    or numerator is zero.
    :param a:  the numerator
    :param b:  the denominator
    :return:   a / b
    """
    if isinstance(a, np.ma.MaskedArray):
        a = a.data
    if isinstance(b, np.ma.MaskedArray):
        b = b.data
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[c == np.inf] = 0
        c[np.isnan(c)] = 0
        c = np.nan_to_num(c)
        c = np.ma.masked_array(c, mask=mask)
    return c


def get_spatial_derivative_cartesian(f, x, y):
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

    dfdx = divide_by_zero(dfdj, dxdj, mask) + divide_by_zero(dfdi, dxdi, mask)
    dfdy = divide_by_zero(dfdi, dydi, mask) + divide_by_zero(dfdj, dydj, mask)

    return dfdx, dfdy


def get_spatial_derivative_cylindrical(f, r, t):
    """
    Returns a derivative of f with respect to r and theta for discreet nd arrays of points.
    :param f:   2d array of some function f values
    :param r:   2d array of radius coordinates (radius)
    :param t:   2d array of tangential coordinates (theta)
    """

    mask = np.zeros(f.shape)
    if isinstance(f, np.ma.MaskedArray):
        mask = f.mask

    # derivatives in terms of index locations i,j
    dfdi, dfdj = np.gradient(f)
    drdi, drdj = np.gradient(x)
    dtdi, dtdj = np.gradient(y)

    dfdr = divide_by_zero(dfdj, drdj, mask) + divide_by_zero(dfdi, drdi, mask)
    dfdt = divide_by_zero(dfdj, dtdj, mask) + divide_by_zero(dfdi, dtdi, mask)

    return dfdr, dfdt


def get_spatial_derivative(f, a, b, cylindrical=False):

    if cylindrical:
        return get_spatial_derivative_cylindrical(f, a, b)
    else:
        return get_spatial_derivative_cartesian(f, a, b)


# testing area
if __name__ == "__main__":

    xl = np.linspace(1, 6, 6)
    yl = np.linspace(1, 6, 6)
    x, y = np.meshgrid(xl, yl)
    r = ((x + 10.5) ** 2 + (y + 10.5) ** 2) ** 0.5
    t = np.arctan(y / x)

    f1 = 3 * x + 4 * y
    dfdx, dfdy = get_spatial_derivative(f1, x, y)
    #print dfdx
    #print dfdy

    f2 = 2*t
    dfdr, dfdt = get_spatial_derivative(f2, r, t, cylindrical=True)
    print r
    print t
    print dfdr
    print dfdt
