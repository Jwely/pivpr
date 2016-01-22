__author__ = "Jwely"

import os
from py.piv import construct_experiments
from py.config import *


def process_experiments(ids=None, min_points=DEFAULT_MIN_POINTS,
                        include_dynamic=False, force_recalc=False):

    exps = construct_experiments(EXPERIMENT_TABLE_PATH, DATA_FULL_DIR,
                                 ids, min_points, include_dynamic, force_recalc)


if __name__ == "__main__":
    process_experiments(include_dynamic=True, force_recalc=False)