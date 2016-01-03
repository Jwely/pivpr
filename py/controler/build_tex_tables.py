__author__ = 'Jwely'

from py.utils import csv_to_tex
import os


def main():
    """
    Execution of this file will build .tex files from .csv files in in the `tables`
    directory. The function calls have been simplified as much as possible to allow
    nicely formatted tables that meet a consistent style.
    """

    table_dir = r"..\..\texdocs\tables"

    csv_to_tex(os.path.join(table_dir, "test_matrix_table.csv"),
               justification="|ccc||ccc||ccc|",
               caption="Experimental conditions for all 70 experiments",
               horizontal_line_rows=[2])

    csv_to_tex(os.path.join(table_dir, "piv_upsampling_displacement.csv"),
               justification="|ccc|",
               caption="Pixel displacements by up sampling order.",
               horizontal_line_rows=[1])

    # the measurements for each station
    for station in [1, 2, 3, 4, 5, 6, 7]:
        csv_to_tex(os.path.join(table_dir, "station_{0}_measurements.csv".format(station)),
                   justification="|ccccccccccc|",
                   caption="Experimental measurements for station {0}".format(station),
                   horizontal_line_rows=[2])

    # v3d format example
    csv_to_tex(os.path.join(table_dir, "v3d_row_example.csv"),
               justification="|cccccccc|",
               caption="Example rows from v3d files with raw 3d vector data",
               horizontal_line_rows=[1])


if __name__ == "__main__":
    main()
