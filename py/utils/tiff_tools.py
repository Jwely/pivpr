__author__ = 'Jwely'

import os
from libtiff import TIFF


def scale_array(array, new_min, new_max, type="linear"):
    """
    simple numeric scaling function to adjust the dynamic range
    of some input set.

    :param array:   array
    :param new_min: new maximum value for input array
    :param new_max: new miniimum value for input array
    :param type: presently only supports linear scaling.
    :return:
    """

    oldtype = array.dtype
    array = array.astype("float")

    if type == "linear":
        old_min = array.min()
        old_max = array.max()

        output = (array + (new_min - old_min)) * (new_max / old_max)
        return output.astype(oldtype)


def load_image_as_array(imagepath):
    """
    loads an image as an array
    """

    tiff = TIFF.open(imagepath, 'r')
    array = tiff.read_image()
    tiff.close()
    return array


def save_array_as_dtype(array, dtype, outpath):
    """
    Saves an array with a specific datatype.

    :param array:   numpy array or array like
    :param dtype:   dtype such as "float32" or "uint16"
    :param outpath: filepath to store tiff
    :return:
    """

    outdir = os.path.dirname(outpath)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # do some integer scaling
    if dtype == "uint16":
        array = scale_array(array, 0, 2 ** 15)
    elif dtype == "uint8":
        array = scale_array(array, 0, 2 ** 7)

    array = array.astype(dtype)
    tiff = TIFF.open(outpath, 'w+')
    tiff.write_image(array)
    tiff.close()