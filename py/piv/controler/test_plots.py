__author__ = 'Jwely'

from py import constructor
from matplotlib import cm


def test_plots(experiment_id):
    """
    This is a scratch board is function for generating plots to your hearts content. mostly just
    used for testing and development
    :param experiment_ids:
    :return:
    """


    exp = constructor.experiments(experiment_table_path="../constructor/dat/experiment_table.csv",
                                  experiment_directory_path="../../data_full",
                                  ids=[experiment_id],
                                  min_points=20,
                                  include_dynamic=True,
                                  force_recalc=False)
    av = exp[0].axial_vortex
    # av.stream_plot()
    # av.quiver_plot()
    # for component in ['U','V','W','u','v','uu','vv','ww','ctke']:
    # av.contour_plot(component)
    # av.contour_plot('t_meshd')
    # av.contour_plot('hv_meshd')
    # av.scatter_plot('r_mesh', 'ctke', 'num', cmap=cm.jet,
    # x_label="Distance to vortex core (mm)",
    # c_label="angle from right horizontal")

    # some ctke plots
    components = ['ctke', 'rt']
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
    test_plots(55)

    import time