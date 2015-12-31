__author__ = 'Jwely'

from py import constructor


def process_all_experiments(force_recalc=True):
    """ Quick one call function to process all datasets, does not return outputs!"""
    constructor.experiments(experiment_table_path="../constructor/dat/experiment_table.csv",
                            experiment_directory_path="../../data_full",
                            ids=None,
                            min_points=20,
                            include_dynamic=True,
                            force_recalc=force_recalc)


if __name__ == "__main__":
    process_all_experiments(force_recalc=True)
