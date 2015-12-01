__author__ = 'Jwely'

from VecFieldCartesian import VecFieldCartesian
from Timer import Timer

import matplotlib.pyplot as plt
import numpy as np
import os
import cPickle


class MeanVecFieldCartesian:

    def __init__(self, v3d_paths, name_tag):
        """
        :param v3d_paths:       list of filepaths to v3d files
        :param name_tag:        unique string name tag for this data set
        :return:
        """

        self.name_tag = name_tag
        self.v3d_paths = v3d_paths
        self.VecFields = []

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
                           'M': None,       # velocity magnitude (all three components)

                           'u': None,       # fluctuation in U
                           'v': None,       # fluctuation in V
                           'w': None,       # fluctuation in W

                           'uu': None,      # turbulent energy in u
                           'vv': None,      # turbulent energy in v
                           'ww': None,      # turbulent energy in w
                           'cte': None,     # total cartesian turbulent energy

                           'uv': None,      # reynolds stress in u/v
                           'uw': None,      # reynolds stress in u/w
                           'vw': None,      # reynolds stress in v/w
                           'crs': None,     # total cartesian reynolds stress

                           'num': None}     # number of good data points making up stats for other values


        # Build up the data set
        self._ingest_paths(v3d_paths)       # creates VecFields from input filepath list
        self._average_cartesian()           # fill in the velocity matrix


    def __getitem__(self, item):
        """ allows components of the vel_matrix to be accessed more simply through instance[key] """
        return self.vel_matrix[item]


    def _ingest_paths(self, filepath_list):
        """ Creates VecField3d objects from each filepath in the list"""
        t = Timer()

        # grab the first one and inherit its dimensional information
        first_vf = VecFieldCartesian(filepath_list.pop(0))
        self.x_set = first_vf.x_set
        self.y_set = first_vf.y_set
        self.dims = first_vf.dims
        self.meshgrid = first_vf.meshgrid
        self.VecFields.append(first_vf)

        for path in filepath_list:
            next_vf = VecFieldCartesian(path)
            assert next_vf.dims == first_vf.dims, "Inconsistent dimensions detected!"
            self.VecFields.append(next_vf)

        print("loaded {0} files in {1} s".format(len(filepath_list), t.finish()))


    def _average_cartesian(self):
        """
        Populates the  U, V, W, u, v, w, uu, vv, ww, uv, uw, and vw  values of the velocity matrix.
        """

        depth = len(self.VecFields)
        u_set = np.ma.zeros(self.dims + tuple([depth]))
        v_set = np.ma.zeros(self.dims + tuple([depth]))
        w_set = np.ma.zeros(self.dims + tuple([depth]))

        print("Taking statistics...")
        for i, vf in enumerate(self.VecFields):
            u_set[:, :, i] = vf.vel_matrix['U']
            v_set[:, :, i] = vf.vel_matrix['V']
            w_set[:, :, i] = vf.vel_matrix['W']

        # populate the velocity matrix
        self.vel_matrix['U'] = np.ma.mean(u_set, axis=2)
        self.vel_matrix['V'] = np.ma.mean(v_set, axis=2)
        self.vel_matrix['W'] = np.ma.mean(w_set, axis=2)
        self.vel_matrix['M'] = (self.vel_matrix['U'] ** 2 +
                                self.vel_matrix['V'] ** 2 +
                                self.vel_matrix['W'] ** 2) ** 0.5

        self.vel_matrix['u'] = np.ma.std(u_set, axis=2)
        self.vel_matrix['v'] = np.ma.std(v_set, axis=2)
        self.vel_matrix['w'] = np.ma.std(w_set, axis=2)

        self.vel_matrix['uu'] = self.vel_matrix['u'] * self.vel_matrix['u']
        self.vel_matrix['vv'] = self.vel_matrix['v'] * self.vel_matrix['v']
        self.vel_matrix['ww'] = self.vel_matrix['w'] * self.vel_matrix['w']
        self.vel_matrix['cte'] = (self.vel_matrix['uu'] + self.vel_matrix['vv'] + self.vel_matrix['ww'])

        self.vel_matrix['uv'] = self.vel_matrix['u'] * self.vel_matrix['v']
        self.vel_matrix['uw'] = self.vel_matrix['u'] * self.vel_matrix['w']
        self.vel_matrix['vw'] = self.vel_matrix['v'] * self.vel_matrix['w']
        self.vel_matrix['crs'] = (self.vel_matrix['uv'] + self.vel_matrix['uw'] + self.vel_matrix['vw'])

        self.vel_matrix['num'] = u_set.count(axis=2)


    def to_pickle(self, pickle_path, reduce_memory=False):
        """ dumps the contents of this object to a pickle """

        # delete the constituent objects with hundreds of additional matrices to reduce pkl size
        if reduce_memory:
            del self.VecFields
            self.VecFields = None

        # create the directory and write the pkl file.
        if not os.path.exists(os.path.dirname(pickle_path)):
            os.mkdir(os.path.dirname(pickle_path))

        with open(pickle_path, 'wb+') as f:
            cPickle.dump(self, f)
            print("Saved to {0}".format(pickle_path))


    @staticmethod
    def from_pickle(pickle_path):
        """ loads previous saved state from a .pkl file and returns a MeanVecFieldCartesian instance """

        with open(pickle_path, 'rb') as f:
            new_instance = cPickle.load(f)

        print("loaded pkl from {0}".format(pickle_path))
        return new_instance


    def show_heatmap(self, component):
        """ prints a quick simple heads up heatmap of input component of the vel_matrix attribute"""
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

    mvf = MeanVecFieldCartesian(paths, "test")
    #mvf.show_heatmap('num')
    print mvf['U']

