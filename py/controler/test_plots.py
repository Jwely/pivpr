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
    #av.quiver_plot()
    #for component in ['U','V','W','u','v','uu','vv','ww','ctke']:
    #av.contour_plot(component)
    #av.contour_plot('t_meshd')
    #av.contour_plot('hv_meshd')
    #av.scatter_plot('r_mesh', 'ctke', 'num', cmap=cm.jet,
    #x_label="Distance to vortex core (mm)",
    #c_label="angle from right horizontal")

    # some ctke plots
    components = ['ctke']
    r_ranges = [('0r', '0.6r'), ('0.6r', '1.4r'), ('1.4r', '2r'), ('2r', '3r')]
    t_ranges = [(0, 90), (0, 45), (45, 90), (0, 30), (30, 60), (60, 90)]
    for component in components:
        for r_range in r_ranges:
            for t_range in t_ranges:
                outname = "exp{0}_dyn_{1}_{2}_{3}".format(experiment_id, component, "-".join(r_range), "-".join(t_range))
                outpath = "../test_plots/{name}.png"
                av.plot_dynamic(component,
                                title=outname,
                                y_label="Turbulent Kinetic Energy",
                                r_range=r_range,
                                t_range=r_range,
                                symmetric=True)


if __name__ == "__main__":
    test_plots(55)

    import time