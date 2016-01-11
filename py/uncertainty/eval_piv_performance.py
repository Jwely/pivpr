__author__ = 'Jwely'

from py.piv import VecFieldCartesian
import numpy as np
from matplotlib import pyplot as plt
import json


def eval_piv_performance(v3d_filepath, params_filepath):
    """
    Evaluates the performance of the PIV system by taking processed vector files from
    images generated with an ArtificialPIV instance and comparing the results against
    the known velocity vectors used to generate the image.
    :param v3d_filepath:
    :return:
    """

    piv_result = VecFieldCartesian(v3d_filepath)
    piv_params = json.loads(open(params_filepath).read())

    # evaluate the u velocity component
    u = piv_result['U'].flatten()
    plt.hist(u[~u.mask], bins=20, align="left", color="darksalmon", normed=True)
    plt.axvline(piv_params['x_vel'], color="navy", linestyle="dashed", linewidth=2)
    plt.show()






if __name__ == "__main__":
    eval_piv_performance("artificial_images/loc-1-high_x/loc-1-high_x.v3d",
                         "artificial_images/loc-1-high_x/params.json",)