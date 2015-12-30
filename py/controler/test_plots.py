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
                                      force_recalc=True)
        av = exp[0].axial_vortex
        #av.stream_plot()
        #av.quiver_plot()
        for component in ['num','r','t','w','rt','rw','tw','ctke']:
            av.contour_plot(component)
        #av.contour_plot('t_meshd')
        #av.contour_plot('hv_meshd')
        #av.scatter_plot('r_mesh', 'ctke', 'num', cmap=cm.jet,
                       #x_label="Distance to vortex core (mm)",
                       #c_label="angle from right horizontal")


if __name__ == "__main__":
    test_plots([55])

    import time
    time.sleep(100)