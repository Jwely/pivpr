__author__ = 'Jwely'

import numpy as np


def masked_rms(ndarray, axis, mask):
    """ quick little function to shorten some often repeated function calls """
    return np.ma.masked_array(np.ma.sqrt(np.ma.mean(ndarray ** 2, axis=axis)), mask=mask)


