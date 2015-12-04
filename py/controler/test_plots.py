__author__ = 'Jwely'

from py import constructor
from matplotlib import cm


def test_plots(experiment_id):
    """
    This is a stratch board is function for generating plots to your hearts content
    :param experiment_id:
    :return:
    """

    exp = constructor.experiments(experiment_table_path="../constructor/dat/experiment_table.csv",
                                  experiment_directory_path="../../data_full",
                                  ids=[experiment_id],
                                  min_points=20,
                                  force_recalc=False)
    av = exp[0].axial_vortex
    #av.stream_plot()
    av.contour_plot('num')
    av.contour_plot('ctke')
    av.contour_plot(['R', 'T', 'W'])
    av.contour_plot(['rr', 'tt', 'ww'])
    #av.contour_plot(['rt', 'rw', 'tw'])
    av.scatter_plot('r_mesh', 'T', 'ctke', cmap=cm.jet,
                   title="Tangential Velocity and Turbulent Kinetic Energy",
                   x_label="Distance to vortex core (mm)",
                   y_label="Tangential velocity (m/s)",
                   c_label="Turbulent kinetic energy (TKE)")


if __name__ == "__main__":
    test_plots(1)