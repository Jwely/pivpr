__author__ = 'Jwely'

import pandas as pd
import os
from py.piv.Experiment import Experiment
from py.config import *
from construct_axial_vortex import construct_axial_vortex


def construct_experiments(experiment_table_path, experiment_directory_path, ids=None,
                          min_points=20, include_dynamic=True, force_recalc=False):
    """
    Constructs an Experiment instance with all useful attributes of the experiment. Some of these
    attributes are read from ancillary data in the `dat` folder.

    :param experiment_table_path:       filepath to `experiment_table.csv` usually in `dat` folder
    :param experiment_directory_path:   path to directory with all exp data `data_full`
    :param ids:                         list of run ID numbers to use
    :param min_points:                  when building the experiments AxialVortex data, cells
                                        with fewer good samples than min_points will be masked and
                                        treated as NoData. A higher min_point requirement reduces
                                        the overall size of the data set, but improves quality
    :param force_recalc:                forces re-computation from raw data of all derived values.
                                        if left False, data may be loaded from a previous binary file
                                        instead of re-computed.
    :return experiments:                a list full of Experiment instances, without dynamic data.
                                        probably pretty memory intensive.
    """

    if isinstance(ids, int):
        ids = [ids]

    if not os.path.exists(experiment_table_path):
        raise Exception("Path does not exist at {0}".format(experiment_table_path))
    dataframe = pd.read_csv(experiment_table_path)
    experiments_list = []

    if ids is None:
        ids = range(0, len(dataframe) + 1)

    for i, row in dataframe.iterrows():
        if row['experiment_id'] in ids:
            # build the experiment object
            kwargs = row.to_dict()
            exp_dir = os.path.join(experiment_directory_path, str(row['experiment_id']))
            exp = Experiment(**kwargs)

            # now build up the vortex associated with it, and add it to the experiment
            name_tag = "ID-{0}_Z-{1}_Vfs-{2}".format(row['experiment_id'],
                                                     row['z_location'],
                                                     row['v_fs_mean'])

            # construct an axial vortex
            av = construct_axial_vortex(v3d_dir=exp_dir,
                                        pkl_dir=PIV_PICKLE_DIR,
                                        name_tag=name_tag,
                                        include_dynamic=include_dynamic,
                                        velocity_fs=row['v_fs_mean'],
                                        eta_p=exp.eta_p,
                                        z_location=int(row['z_location'] * 25.4 + 0.5),     # rounded z location in mm
                                        force_recalc=force_recalc,
                                        min_points=min_points)
            exp.ingest_axial_vortex(av)
            exp.to_json(os.path.join(PIV_PICKLE_DIR, "{0}.json".format(name_tag)))
            experiments_list.append(exp)

    return experiments_list


if __name__ == "__main__":
    construct_experiments(EXPERIMENT_TABLE_PATH, DATA_FULL_DIR)




