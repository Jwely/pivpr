__author__ = 'Jwely'

from py import constructor
import os


class FinalFigure:
    """
    The template class for generating final figures from experimental datasets.
    Anything that extends this class is producing final figures which are included in
    the final thesis document which can be found at [github.com/Jwely/pivpr_docs]

    All extensions should have a "generate"  method, and require no arguments upon initiation,
    this allows simple systematic generation of figures.
    """

    def __init__(self, experiment_ids, sub_dir):

        if not isinstance(experiment_ids, list):
            experiment_ids = [experiment_ids]

        self.experiment_number = experiment_ids
        self.fig_dir = os.path.abspath("../final_figures/{0}".format(sub_dir))
        self.experiments = constructor.experiments(
            experiment_table_path="../constructor/dat/experiment_table.csv",
            experiment_directory_path="../../data_full",
            ids=experiment_ids,
            min_points=20,
            force_recalc=False)
