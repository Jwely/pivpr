__author__ = 'Jwely'

from py.tex import csv_to_tex
from py.config import *
import os


def build_tex_tables():
    """
    Execution of this file will build .tex files from .csv files in in the `tables`
    directory. The function calls have been simplified as much as possible to allow
    nicely formatted tables that meet a consistent style.
    """

    csv_to_tex(os.path.join(TEX_TABLE_DIR, "test_matrix_table.csv"),
               justification="|ccc||ccc||ccc|",
               caption="Experimental conditions for all 70 experiments",
               horizontal_line_rows=[2])

    csv_to_tex(os.path.join(TEX_TABLE_DIR, "piv_upsampling_displacement.csv"),
               caption="Pixel displacements by up sampling order.",
               horizontal_line_rows=[1])

    # the measurements for each station
    for station in [1, 2, 3, 4, 5, 6, 7]:
        csv_to_tex(os.path.join(TEX_TABLE_DIR, "station_{0}_measurements.csv".format(station)),
                   caption="Experimental conditions for tests at station {0}".format(station),
                   horizontal_line_rows=[2])

    # v3d format example
    csv_to_tex(os.path.join(TEX_TABLE_DIR, "v3d_row_example.csv"),
               caption="Example rows from v3d files with raw 3d vector data",
               horizontal_line_rows=[1])

    # experiment results summary table
    csv_to_tex(os.path.join(TEX_TABLE_DIR, "experiment_results_summary.csv"),
               caption="Summary of Experimental Results",
               horizontal_line_rows=[1])

    # table of uncertainty simulation scenarios
    csv_to_tex(os.path.join(TEX_TABLE_DIR, "uncertainty_sim_table.csv"),
               caption="Velocity conditions of Monte Carlo image generation for uncertainty analysis.",
               horizontal_line_rows=[2])


if __name__ == "__main__":
    build_tex_tables()
