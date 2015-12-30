__author__ = 'Jwely'

from py import constructor
from matplotlib import cm


def test_plots(experiment_ids):
    """
    This is a scratch board is function for generating plots to your hearts content
    :param experiment_ids:
    :return:
    """

    for experiment_id in experiment_ids:
        exp = constructor.experiments(experiment_table_path="../constructor/dat/experiment_table.csv",
                                      experiment_directory_path="../../data_full",
                                      ids=[experiment_id],
                                      min_points=20,
                                      force_recalc=False)
        av = exp[0].axial_vortex
        #av.stream_plot()
        #av.quiver_plot()
        #av.contour_plot('P')
        av.contour_plot('t_meshd')
        av.contour_plot('hv_meshd')
        #v.contour_plot(['R', 'T', 'W'])
        #av.contour_plot(['rr', 'tt', 'ww'])
        #av.contour_plot(['rt', 'rw', 'tw'])
        av.scatter_plot('r_mesh', 'ctke', 'hv_meshd', cmap=cm.jet,
                       x_label="Distance to vortex core (mm)",
                       c_label="angle from right horizontal")


if __name__ == "__main__":
    test_plots([55])

    import time
    time.sleep(100)