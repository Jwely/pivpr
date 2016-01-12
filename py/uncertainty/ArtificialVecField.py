__author__ = 'Jwely'

import json
import numpy as np
from matplotlib import pyplot as plt
from py.piv import VecFieldCartesian



class ArtificialVecField(VecFieldCartesian):
    """
    Class for evaluating performance of PIV to resolve images synthesized
    by the ArtificialPIV class with precisely known particle movements. It
    extends the VecFieldCartesian class, and adds just a couple methods
    for gathering error and uncertainty information.
    """

    def __init__(self, v3d_path, params_filepath):
        """
        :param v3d_path: filepath to v3d file
        :param params_filepath: filepath to params.json file
        :return:
        """
        VecFieldCartesian.__init__(self, v3d_path)
        self.piv_params = json.loads(open(params_filepath).read())


    def plot_histogram(self, component, title=None, outpath=None):
        """
        Creates a histogram of the vectors within the artificial vector field
        and compares it against the parameters used to generate the vector field.
        Thus, bias and random error specific to the parameters used to generate
        the artificial PIV images can be obtained.

        :param component:   either 'U', 'V', or 'W'
        :param title:       custom title for the plot
        :param outpath:     filepath to save the figure to.
        :return:
        """

        figsize = (10, 5)    # hardcoded for now

        data = self[component].flatten()

        results = {"top": np.nanpercentile(data.filled(np.nan), 99.85),
                   "bot": np.nanpercentile(data.filled(np.nan), 0.15),
                   "mean": np.mean(data),
                   "sim": self.piv_params[component.lower()]}

        results.update({"top_er": 100 * (results['top'] - results['sim']) / results['sim'],
                        "bot_er": 100 * (results['bot'] - results['sim']) / results['sim'],
                        "mean_er": 100 * (results['mean'] - results['sim']) / results['sim']})

        # create the histogram
        fig = plt.figure(figsize=figsize, dpi=120, facecolor='w', edgecolor='k')
        plt.hist(data[~data.mask], bins=100, align="left", color="grey", alpha=0.4)

        # actual velocity line
        plt.axvline(results["sim"], color="navy",
                    linestyle="-", linewidth=1, label="Simulated")

        # percentile and average lines
        plt.axvline(results["mean"], color="navy", linestyle="--", linewidth=1, label="Measured Mean")
        plt.axvline(results["top"], color="navy", linestyle=":", linewidth=1, label="99.7% Confidence")
        plt.axvline(results["bot"], color="navy", linestyle=":", linewidth=1)

        plt.xlabel("${0}$  $(m/s)$".format(component.lower()))
        plt.ylabel("Number of Occurrences")
        plt.legend(loc=2)

        if title is not None:
            plt.title(title)

        if outpath:
            plt.savefig(outpath)
            print("saved figure to {0}".format(outpath))
            plt.close()
        else:
            plt.show(fig)

        return results


if __name__ == "__main__":
    avf = ArtificialVecField("artificial_images/loc-1-high_x/loc-1-high_x.v3d",
                             "artificial_images/loc-1-high_x/params.json",)
    #avf.plot_histogram('U')
    #avf.plot_histogram('V')
    print avf.plot_histogram('W', "Histogram of $w$")