__author__ = 'Jwely'

from py.piv.construct_experiments import construct_experiments
from py.config import *


def process_all_experiments(force_recalc=True):
    """ Quick one call function to process all datasets, does not return outputs!"""
    construct_experiments(experiment_table_path=EXPERIMENT_TABLE_PATH,
                          experiment_directory_path=DATA_FULL_DIR,
                          ids=None,
                          min_points=20,
                          include_dynamic=True,
                          force_recalc=force_recalc)


if __name__ == "__main__":
    process_all_experiments(force_recalc=False)
