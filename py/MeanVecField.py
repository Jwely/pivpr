__author__ = 'Jwely'

from VecField import VecField
from Timer import Timer

import numpy as np


class MeanVecField3d:


    def __init__(self, v3d_paths, name_tag):
        """

        :param v3d_paths:       list of filepaths to v3d files
        :param name_tag:        unique string name tag for this data set
        :return:
        """

        self.name_tag = name_tag
        self.v3d_paths = v3d_paths
        self.VecFields = self._ingest_paths(v3d_paths)

        # dictionary of matrices by key symbol. Capitols are averages, lowercase are fluctuations
        self.velocity_matrix = {'U': None,      # x direction mean velocity
                                'V': None,      # y direction mean velocity
                                'W': None,      # z direction mean velocity
                                'R': None,      # mean radial velocity around vortex core
                                'T': None,      # mean tangential velocity around vortex core

                                'u': None,      # fluctuation in U
                                'v': None,      # fluctuation in V
                                'w': None,      # fluctuation in W
                                'r': None,      # fluctuation in R
                                't': None,      # fluctuation in T

                                'uu': None,     # turbulent energy in u
                                'vv': None,     # turbulent energy in v
                                'uv': None,     # reynolds stress in u/v
                                'uw': None,     # reynolds stress in u/w
                                'vw': None,     # reynolds stress in v/w

                                'rr': None,     # turbulent energy in radial
                                'tt': None,     # turbulent energy in tangential
                                'rt': None,     # reynolds stress in r/t
                                'rw': None,     # reynolds stress in r/w
                                'tw': None,     # reynolds stress in t/w

                                'ww': None,     # turbulent energy in w
                                'tot': None,    # total turbulent energies
                                'num': None}    # number of good data points making up stats for other values

        # all the same data in the matrices above, but flattened into a 1d array and in terms of "distance to core"
        self.velocity_flat = {'U': None,        # x direction mean velocity
                              'V': None,        # y direction mean velocity
                              'W': None,        # z direction mean velocity
                              'R': None,        # radial mean velocity around vortex core
                              'T': None,        # tangential mean velocity around vortex core

                              'u': None,        # U fluctuation
                              'v': None,        # V fluctuation
                              'w': None,        # W fluctuation
                              'r': None,        # R fluctuation
                              't': None,        # T fluctuation

                              'uu': None,       # turbulent energy in u
                              'vv': None,       # turbulent energy in v
                              'uv': None,       # reynolds stress in u/v
                              'uw': None,       # reynolds stress in u/w
                              'vw': None,       # reynolds stress in v/w

                              'rr': None,       # turbulent energy in radial
                              'tt': None,       # turbulent energy in tangential
                              'rt': None,       # reynolds stress in r/t
                              'rw': None,       # reynolds stress in r/w
                              'tw': None,       # reynolds stress in t/w

                              'ww': None,       # turbulent energy in w
                              'tot': None,      # total turbulent energies
                              'num': None}      # number of good data points making up stats for other values


    @staticmethod
    def _ingest_paths(filepath_list):
        """ Creates VecField3d objects from each filepath in the list"""
        t = Timer()
        outlist = []

        # use of "set" ensures no duplicate file paths
        for path in set(filepath_list):
            outlist.append(VecField(path))

        print("loaded {0} files in {1} s".format(len(filepath_list), t.finish()))
        return outlist


    def to_pickle(self):
        pass


    def from_pickle(self):
        pass


    def average_cartesian(self, *args):
        """
        Populates the  U, V, W, u, v, w, uu, vv, ww, uv, uw, and vw  values of the velocity matrix.
        :param args:
        :return:
        """

        depth = len(self.VecFields)

        for vf3d in self.VecFields:





    def find_vortex_core(self):
        pass





if __name__ == "__main__":

    paths = [r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01000.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01001.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01002.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01003.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01004.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01005.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01006.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01007.v3d"]

    mvf = MeanVecField3d(paths, "test")
