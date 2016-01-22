__author__ = "Jwely"

import numpy as np


def get_dydx(x, y):
    """
    simply returns a step by step arithmetic derivative of dy/dx.
    input two array like series of identical length. Function will return
    two arrays which are one item shorter than the inputs. This is because the
    dydx values are sampled at the midpoints of each line segement.

    :param x:           x array
    :param y:           y array
    :return: xdx, dydx
    """

    dydx = np.zeros(len(x) - 1)
    xdx = np.zeros(len(x) - 1)

    i = 0
    while i < len(x) - 1:
        dx = x[i + 1] - x[i]
        dy = y[i + 1] - y[i]
        dydx[i] = dy / dx
        xdx[i] = (x[i] + x[i + 1]) / 2
        i += 1

    return xdx, dydx



