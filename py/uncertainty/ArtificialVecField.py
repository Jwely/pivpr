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
        self.piv_params = json.loads(open(params_filepath).read())
        VecFieldCartesian.__init__(self, v3d_path, velocity_fs=self.piv_params['w'])
        self.error_data = {}

        # correction factor for directional convention in Y direction
        self.vel_matrix['V'] *= -1


    def get_error(self, component, n_measurements):

        # compile some numerical results from the data to return.
        data = self[component].flatten()
        results = {"top": np.nanpercentile(data.filled(np.nan), 97.5),
                   "bot": np.nanpercentile(data.filled(np.nan), 2.5),
                   "mean": np.mean(data),
                   "sim": self.piv_params[component.lower()]}

        # update results to be in terms of bias and precision for a total uncertainty
        num = data.count()
        delta = data - self.piv_params[component.lower()]
        bias = np.mean(delta)
        stdev = np.ma.std(data)
        precision = stdev * 2 / (n_measurements ** 0.5)
        uncertainty = (bias ** 2 + precision ** 2) ** 0.5

        results.update({"top_dif": results['top'] - results['sim'],
                        "bot_dif": results['bot'] - results['sim'],
                        "bias": bias,
                        "stdev": stdev,
                        "precision": precision,
                        "uncertainty": uncertainty})

        self.error_data[component] = results
        return results


    def plot_histogram(self, component, n_measurements=200.0, title=None, outpath=None):
        """
        Creates a histogram of the vectors within the artificial vector field
        and compares it against the parameters used to generate the vector field.
        Thus, bias and random error specific to the parameters used to generate
        the artificial PIV images can be obtained.

        :param component:       either 'U', 'V', or 'W'
        :param n_measurements:  The number of measurements taken to determine undertainty (usually 200)
        :param title:           custom title for the plot
        :param outpath:         filepath to save the figure to.
        :return:
        """

        figsize = (8, 2.5)    # hardcoded for now

        # compile some numerical results from the data to return.
        data = self[component].flatten()
        results = self.get_error(component, n_measurements)

        # create the histogram
        fig = plt.figure(figsize=figsize, dpi=120, facecolor='w', edgecolor='k')
        plt.hist(data[~data.mask], bins=150, align="left", color="grey", alpha=0.4)
        plt.xlim(results['mean'] - 3 * results['stdev'],
                 results['mean'] + 4 * results['stdev'])

        # actual velocity line
        plt.axvline(results["sim"], color="navy",
                    linestyle="-", linewidth=1, label="Simulated")

        # percentile and average lines
        plt.axvline(results["mean"], color="navy", linestyle="--", linewidth=1, label="Measured Mean")
        plt.axvline(results["top"], color="navy", linestyle=":", linewidth=1, label="95% Confidence")
        plt.axvline(results["bot"], color="navy", linestyle=":", linewidth=1)

        plt.xlabel("${0}$  $(m/s)$".format(component.lower()))
        plt.ylabel("Frequency")
        plt.legend(loc=0)
        plt.tight_layout()

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
    avf = ArtificialVecField("artificial_images/Ely_May28th07001.v3d",
                             "artificial_images/Ely_May28th07001.json",)

    print avf.plot_histogram('U')
    print avf.plot_histogram('V')
    print avf.plot_histogram('W', 200, title="Histogram of $w$")