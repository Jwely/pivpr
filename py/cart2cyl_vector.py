__author__ = 'Jwely'

from numpy import math as npm


def cart2cyl_vector(u, v, t_mesh):
    """

    :param u:       set of x component of vector
    :param v:       set of y component of vector
    :param t_mesh:  meshgrid for tangential coordinates (t)
    :return: r, t
    """

    r = u * npm.cos(t_mesh) + v * npm.sin(t_mesh)
    t = v * npm.cos(t_mesh) - u * npm.sin(t_mesh)
    return r, t

