__author__ = "Jwely"

import numpy


def smooth_filt(y, window_len=11, window='hanning'):
    """
    smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    :param y:           the input signal
    :param window_len:  the dimension of the smoothing window; should be an odd integer > 5
    :param window:      the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
                        flat window will produce a moving average smoothing.
    :return:            the smoothed signal

    see also:
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
        numpy.convolve
    """

    if y.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if y.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return y

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window may be 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = numpy.r_[y[window_len - 1:0:-1], y, y[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval("numpy.{0}(window_len)".format(window))

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y[(window_len / 2 - 1):-(window_len / 2 + 1)]
