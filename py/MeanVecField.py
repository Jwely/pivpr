__author__ = 'Jwely'

from VecField import VecField
from Timer import Timer

import matplotlib.pyplot as plt
import numpy as np
import os
import pickle


class MeanVecField:

    def __init__(self, v3d_paths, name_tag, pkl_dir=None):
        """
        :param v3d_paths:       list of filepaths to v3d files
        :param name_tag:        unique string name tag for this data set
        :return:
        """

        self.name_tag = name_tag
        self.v3d_paths = v3d_paths
        self.pkl_dir = pkl_dir
        self.pkl_path = None
        self.pkl_path_small = None
        self.VecFields = []

        # establish pkl filepaths to load pre-existing calculations
        self._get_pkl_names()

        # attributes that will be inherited from the constituent VecField instances
        self.x_set = None
        self.y_set = None
        self.dims = (None, None)
        self.meshgrid = {"x": None,
                         "y": None}

        # dictionary of matrices by key symbol. Capitols are averages, lowercase are fluctuations
        self.vel_matrix = {'U': None,       # x direction mean velocity
                           'V': None,       # y direction mean velocity
                           'W': None,       # z direction mean velocity

                           'u': None,       # fluctuation in U
                           'v': None,       # fluctuation in V
                           'w': None,       # fluctuation in W

                           'uu': None,      # turbulent energy in u
                           'vv': None,      # turbulent energy in v
                           'uv': None,      # reynolds stress in u/v
                           'uw': None,      # reynolds stress in u/w
                           'vw': None,      # reynolds stress in v/w

                           'ww': None,      # turbulent energy in w
                           'tot': None,     # total turbulent energies
                           'num': None}     # number of good data points making up stats for other values


        # Build up the data set
        self._ingest_paths(v3d_paths)       # creates VecFields from input filepath list
        self._average_cartesian()           # fill in the velocity matrix


    def _ingest_paths(self, filepath_list):
        """ Creates VecField3d objects from each filepath in the list"""
        t = Timer()

        # grab the first one and inherit its dimensional information
        first_vf = VecField(filepath_list.pop(0))
        self.x_set = first_vf.x_set
        self.y_set = first_vf.y_set
        self.dims = first_vf.dims
        self.meshgrid = first_vf.meshgrid
        self.VecFields.append(first_vf)

        for path in filepath_list:
            next_vf = VecField(path)
            assert next_vf.dims == first_vf.dims, "Inconsistent dimensions detected!"
            self.VecFields.append(next_vf)

        print("loaded {0} files in {1} s".format(len(filepath_list), t.finish()))


    def _get_pkl_names(self):

        if self.pkl_dir is not None:
            self.pkl_path = os.path.join(self.pkl_dir, "{0}_small.pkl".format(self.name_tag))
            self.pkl_path_small = os.path.join(self.pkl_dir, "{0}_small.pkl".format(self.name_tag))


    def _average_cartesian(self):
        """
        Populates the  U, V, W, u, v, w, uu, vv, ww, uv, uw, and vw  values of the velocity matrix.
        """

        depth = len(self.VecFields)
        Uset = np.ma.zeros(self.dims + tuple([depth]))
        Vset = np.ma.zeros(self.dims + tuple([depth]))
        Wset = np.ma.zeros(self.dims + tuple([depth]))

        print("Taking statistics...")
        for i, vf in enumerate(self.VecFields):
            Uset[:, :, i] = vf.vel_matrix['U']
            Vset[:, :, i] = vf.vel_matrix['V']
            Wset[:, :, i] = vf.vel_matrix['W']

        # populate the velocity matrix
        self.vel_matrix['U'] = np.ma.mean(Uset, axis=2)
        self.vel_matrix['V'] = np.ma.mean(Vset, axis=2)
        self.vel_matrix['W'] = np.ma.mean(Wset, axis=2)

        self.vel_matrix['u'] = np.ma.std(Uset, axis=2)
        self.vel_matrix['v'] = np.ma.std(Vset, axis=2)
        self.vel_matrix['w'] = np.ma.std(Wset, axis=2)

        self.vel_matrix['uu'] = self.vel_matrix['u'] * self.vel_matrix['u']
        self.vel_matrix['vv'] = self.vel_matrix['v'] * self.vel_matrix['v']
        self.vel_matrix['ww'] = self.vel_matrix['w'] * self.vel_matrix['w']
        self.vel_matrix['cte'] = self.vel_matrix['uu'] + self.vel_matrix['vv'] + self.vel_matrix['ww']

        self.vel_matrix['uv'] = self.vel_matrix['u'] * self.vel_matrix['v']
        self.vel_matrix['uw'] = self.vel_matrix['u'] * self.vel_matrix['w']
        self.vel_matrix['vw'] = self.vel_matrix['v'] * self.vel_matrix['w']
        self.vel_matrix['crs'] = self.vel_matrix['uv'] + self.vel_matrix['uw'] + self.vel_matrix['vw']

        self.vel_matrix['num'] = Uset.count(axis=2)


    def to_pickle(self, reduce_memory=False):
        """
        dumps the contents of this object to a pickle
        """

        # delete the constituent objects with hundreds of additional matrices to reduce pkl size
        if reduce_memory:
            del self.VecFields
            self.VecFields = None

        # write the pickle
        if not os.path.exists(self.pkl_dir):
            os.mkdir(self.pkl_dir)


    def from_pickle(self):

        if os.path.exists(self.pkl_path) or os.path.exists(self.pkl_path_small):
            # load the pickle here
            pass


    def show(self, component):
        """ prints a quick simple heads up  heatmap of each of the components """
        fig, ax = plt.subplots()
        heatmap = ax.pcolor(self.vel_matrix[component])
        plt.show()



if __name__ == "__main__":
    paths = [r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01000.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01001.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01002.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01003.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01004.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01005.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01006.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01007.v3d"]

    mvf = MeanVecField(paths, "test")
    mvf.show('num')
