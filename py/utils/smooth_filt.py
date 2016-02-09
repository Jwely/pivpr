__author__ = "Jwely"

from scipy import interpolate, signal
import numpy as np


def smooth_filt(x, y, numpoints, M, std, convolve_mode='same', order=1):
    """
    smooths noisy point data
    :param x:               x series
    :param y:               y series
    :param numpoints:       length of vectors, proportional to resolution
    :param M:               number of points in gaussian resample
    :param std:             standard deviation of points in gaussian resample
    :param convolve_mode:   mode of convolution via signal.convolve
    :return:
    """

    def smooth_filt_order1(x, y, numpoints, M=None, std=None, convolve_mode='same', order=None):
        # create interpolation function
        f = interpolate.interp1d(x, y)

        # make linespace
        xx = np.linspace(np.ma.min(x), np.ma.max(x), numpoints)

        # create new yy based on xx linespace
        yy = f(xx)

        # make a window
        window = signal.gaussian(M, std)

        # convolve the signals
        y_smooth = signal.convolve(yy, window / window.sum(), mode='same')

        return xx, y_smooth

    # execute first order
    xx, y_smooth = smooth_filt_order1(x, y, numpoints, M, std, convolve_mode)

    if order == 1:
        return xx, y_smooth
    else:
        for i in xrange(order - 1):
            xx, y_smooth = smooth_filt_order1(xx, y_smooth, numpoints, M, std, convolve_mode)
        return xx, y_smooth



