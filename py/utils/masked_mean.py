__author__ = 'Jwely'

import numpy as np


def masked_mean(ndarray, axis, mask):
    """ quick little function to shorten some often repeated function calls """
    return np.ma.masked_array(np.ma.mean(ndarray, axis=axis), mask=mask)


