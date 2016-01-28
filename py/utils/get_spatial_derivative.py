__author__ = "Jwely"

import numpy as np


def get_spatial_derivative(f, x, y):
    """
    Returns a derivative of f with respect to x for discreet nd arrays of points.
    :param f:   2d array of some function f values
    :param x:   2d array of x coordinate values
    :param y:   2d array of y coordinate values
    """
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

    mask = np.zeros(f.shape)
    if isinstance(f, np.ma.MaskedArray):
        mask = f.mask

    # derivatives in terms of index locations i,j
    dfdi, dfdj = np.gradient(f)
    dxdi, dxdj = np.gradient(x)
    dydi, dydj = np.gradient(y)

    # NOTE: the x and y meshgrids must be perfect, otherwise this function will not deliver complete
    # results at the moment. This means the meshgrids must only actually vary in one dimension, and
    # be simplifiable to 1d arrays. Rotated meshes will not work until the masking issue is solved.
    dfdx = divide_by_zero(dfdj, dxdj, mask)  # + divide_by_zero(dfdi, dxdi, mask)
    dfdy = divide_by_zero(dfdi, dydi, mask)  # + divide_by_zero(dfdj, dydj, mask)
    return dfdx, dfdy


# testing area
if __name__ == "__main__":

    f = np.array([[1, 2, 3],
                  [2, 3, 4],
                  [3, 4, 5]])

    x = np.array([[10, 20, 30],
                  [10, 20, 30],
                  [10, 20, 30]])

    y = np.array([[2, 2, 2],
                  [4, 4, 4],
                  [8, 8, 8]])

    dfdx, dfdy = get_spatial_derivative(f, x, y)




