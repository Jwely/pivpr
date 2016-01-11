__author__ = 'Jwely'

import json
import numpy as np
from matplotlib import pyplot as plt
from py.piv import VecFieldCartesian



class ArtificialVecField(VecFieldCartesian):

    def __init__(self, v3d_path, params_filepath):
        VecFieldCartesian.__init__(self, v3d_path)

        self.piv_params = json.loads(open(params_filepath).read())


    def plot_histogram(self, component):

        data = self[component].flatten()
        plt.hist(data[~data.mask], bins=30, align="left", color="darksalmon", normed=True)
        plt.axvline(self.piv_params[component.lower()], color="navy", linestyle="dashed", linewidth=2)
        plt.show()


if __name__ == "__main__":
    avf = ArtificialVecField("artificial_images/loc-1-high_x/loc-1-high_x.v3d",
                             "artificial_images/loc-1-high_x/params.json",)
    avf.plot_histogram('U')