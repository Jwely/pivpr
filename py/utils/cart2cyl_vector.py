__author__ = 'Jwely'

import numpy as np


def cart2cyl_vector(u, v, t_mesh):
    """
    :param u:       set of x component of vector
    :param v:       set of y component of vector
    :param t_mesh:  meshgrid for tangential coordinates (t)
    :return: r, t
    """

    r = u * np.cos(t_mesh) + v * np.sin(t_mesh)
    t = u * np.sin(t_mesh) - v * np.cos(t_mesh)
    return r, t

