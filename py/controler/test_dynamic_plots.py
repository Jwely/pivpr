__author__ = 'Jwely'

from py.piv.construct_experiments import construct_experiments
from matplotlib import cm
from py.config import *


def test_dynamic_plots(experiment_id):
    """
    This is a scratch board is function for generating plots to your hearts content. mostly just
    used for testing and development
    :param experiment_ids:
    :return:
    """

    exp = construct_experiments(experiment_table_path=EXPERIMENT_TABLE_PATH,
                                experiment_directory_path=DATA_FULL_DIR,
                                ids=[experiment_id],
                                min_points=20,
                                include_dynamic=True,
                                force_recalc=False)
    av = exp[0].axial_vortex

    # some plots!
    components = ['ctke', 'rt', 'w']
    r_ranges = [('0.0r', '0.6r'), ('0.6r', '1.4r'), ('1.4r', '2.0r'), ('2.0r', '3.0r')]
    t_ranges = [(0, 90), (0, 30), (30, 60), (60, 90)]

    for component in components:
        for r_range in r_ranges:
            t_range = (0, 90)
            outname = "exp{0}_dyn_{1}_{2}_{3}".format(experiment_id, component,
                                                      "-".join(r_range), "-".join(map(str, t_range)))
            outpath = "../test_plots/{name}.png".format(name=outname)
            av.dynamic_plot(component,
                            title="",
                            r_range=r_range,
                            t_range=t_range,
                            symmetric=True,
                            outpath=outpath
            )

        for t_range in t_ranges:
            r_range = ('0r', '2r')
            outname = "exp{0}_dyn_{1}_{2}_{3}".format(experiment_id, component,
                                                      "-".join(r_range), "-".join(map(str, t_range)))
            outpath = "../test_plots/{name}.png".format(name=outname)
            av.dynamic_plot(component,
                            title="",
                            r_range=r_range,
                            t_range=t_range,
                            symmetric=True,
                            outpath=outpath
            )


if __name__ == "__main__":
    test_dynamic_plots(55)

    import time