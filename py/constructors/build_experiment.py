__author__ = 'Jwely'

import pandas as pd
import os
from py.managers.Experiment import Experiment
from build_axial_vortex import build_axial_vortex


def build_experiments(experiment_table_path, experiment_directory_path, ids=None):
    """
    Constructs an Experiment instance with all useful attributes of the experiment. Some of these
    attributes are read from ancillary data in the `dat` folder.

    :param experiment_table_path:
    :param experiment_directory_path:
    :param ids:                         list of run ID numbers to use
    :return:
    """

    dataframe = pd.read_csv(experiment_table_path)
    experiments = []

    if ids is None:
        ids = range(0, len(dataframe) + 1)

    for i, row in dataframe.iterrows():
        if row['experiment_id'] in ids:
            # build the experiment object
            kwargs = row.to_dict()
            exp_dir = os.path.join(experiment_directory_path, "1")
            exp = Experiment(experiment_dir=exp_dir, **kwargs)

            # now build up the vortex associated with it, and add it to the experiment
            name_tag = "ID-{0}_Z-{1}_Vfs-{2}".format(row['experiment_id'],
                                                     row['z_location'],
                                                     row['v_fs_mean'])

            av = build_axial_vortex(exp_dir, "../pickles",
                                    name_tag=name_tag,
                                    include_dynamic=False,
                                    velocity_fs=row['v_fs_mean'],
                                    force_recalc=False)
            exp.ingest_axial_vortex(av)
            experiments.append(exp)

    return experiments


if __name__ == "__main__":

    build_experiments("dat/experiment_table.csv", "../../data_full")




