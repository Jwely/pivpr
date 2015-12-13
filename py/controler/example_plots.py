__author__ = 'Jwely'

from py import constructor
from matplotlib import cm


def example_plots(experiment_id):
    """
    Function used for making test case plots
    """

    exp = constructor.experiments(experiment_table_path="../constructor/dat/experiment_table.csv",
                                  experiment_directory_path="../../data_full",
                                  ids=[experiment_id],
                                  min_points=20,
                                  force_recalc=False)
    av = exp[0].axial_vortex
    av.stream_plot()
    for component in ['num', 'T', 'R', 'W', 'rr', 'tt', 'ww', 'ctke', 'rt', 'rw', 'tw']:
        av.contour_plot(component)

    av.scatter_plot('r_mesh', 'T', 'ctke', cmap=cm.jet,
                    title="Tangential Velocity and Turbulent Kinetic Energy",
                    x_label="Distance to vortex core (mm)",
                    y_label="Tangential velocity (m/s)",
                    c_label="Turbulent kinetic energy (TKE)")


if __name__ == "__main__":
    example_plots(59)