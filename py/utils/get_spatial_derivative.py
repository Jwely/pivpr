__author__ = "Jwely"

import numpy


def get_spatial_derivative(f, x):
    """
    Returns a derivative of f with respect to x for discreet nd arrays of points.
    :param f:   nd_array of function values
    :param x:   nd_array of spatial coordinate values

    Example:
        given a matrix of values (f), and meshgrids of coordinates in real
        space (x, y). One may compute partial derivative df/dx with

            get_spatial_derivative(f, x)

        and the partial df/dy with

            get_spatial derivative(f, y)

    :return: df/dx
    """
    df = numpy.gradient(f)
    dx = numpy.gradient(x)
    dfdx = numpy.divide(df, dx)
    return dfdx


