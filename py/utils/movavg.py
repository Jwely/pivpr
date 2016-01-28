__author__ = "Jwely"

import numpy as np


def movavg(x, y, window_size, fixed_density=False):
    """
    returns a moving average of y along x axis within the window size. x and y
    must be one dimensional array-like of identical size.

    :param x:               x series
    :param y:               y series
    :param window_size:     Window width in units of x axis
    :param irregular:       Set to True if x points are not at regular intervals
    :return:
    """

    # ensure both x and y are ordered ascendingly
    ordering = np.argsort(x)
    x = x[ordering]
    y = y[ordering]
    window_size = float(window_size)

    # set up yavg initially as y
    yavg = y

    if fixed_density is False:
        for i, yi in enumerate(yavg):
            irange = np.greater(x, x[i] - window_size / 2) * np.less(x, x[i] + window_size / 2)
            yavg[i] = np.mean(y[irange])
        return x, yavg
    else:
        xfd = np.array(range(0, int(abs(x.data[-1] - x.data[0]) / window_size))) * window_size
        yavgfd = np.zeros(xfd.shape)
        for i, yavgfdi in enumerate(yavgfd[:-1]):
            irange = np.greater(x, xfd[i]) * np.less(x, xfd[i + 1])
            yavgfd[i] = np.mean(y[irange])
        mask = np.isnan(yavgfd)
        return np.ma.masked_array(xfd, mask), np.ma.masked_array(yavgfd, mask)









