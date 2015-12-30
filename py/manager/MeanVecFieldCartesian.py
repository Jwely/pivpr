__author__ = 'Jwely'

import os
import numpy as np
from py.manager.VecFieldCartesian import VecFieldCartesian
from py.utils import Timer
from py.config import *


def masked_rms(ndarray, axis, mask):
    """ quick little function to shorten some often repeated function calls """
    return np.ma.masked_array(np.ma.sqrt(np.ma.mean(ndarray ** 2, axis=axis)), mask=mask)



class MeanVecFieldCartesian:

    def __init__(self, name_tag=None, v3d_paths=None, min_points=20, velocity_fs=None):
        """
        :param v3d_paths:       list of filepaths to v3d files
        :param name_tag:        unique string name tag for this data set
        :param min_points:      minimum number of values for given grid location for statistics to be
                                considered valid. The more points the better, especially for turbulence values.
        :param velocity_fs:     free stream velocity. not really used yet.
        :return:
        """

        self.name_tag = name_tag
        self.v3d_paths = v3d_paths
        self.min_points = min_points
        self.velocity_fs = velocity_fs
        self.constituent_vel_matrix_list = []    # empty list for constituent VecFieldCartesian instances

        # attributes that will be inherited from the constituent VecField instances
        self.x_set = None
        self.y_set = None
        self.dims = (None, None, None)
        self.meshgrid = {"x_mesh": None,
                         "y_mesh": None}

        # dictionary of matrices by key symbol. Capitols are averages, lowercase are fluctuations
        self.vel_matrix = {'U': None,       # x direction mean velocity
                           'V': None,       # y direction mean velocity
                           'W': None,       # z direction mean velocity
                           'P': None,       # in-plane velocities (U and V components, not W)
                           'M': None,       # velocity magnitude (all three components)

                           'u': None,       # mean absolute fluctuation in U
                           'v': None,       # mean absolute fluctuation in V
                           'w': None,       # mean absolute fluctuation in W

                           'uu': None,      # turbulent energy in u (u' * u') bar
                           'vv': None,      # turbulent energy in v (v' * v') bar
                           'ww': None,      # turbulent energy in w (w' * w') bar
                           'ctke': None,    # cartesian turbulent kinetic energy (TKE)

                           'uv': None,      # reynolds stress in u/v
                           'uw': None,      # reynolds stress in u/w
                           'vw': None,      # reynolds stress in v/w

                           'num': None}     # number of good data points making up stats for other values

        # the same as vel_matrix, except every every time slice is included for an n+1 dimensional matrix
        self.dynamic_set = {'U': None,       # x direction mean velocity
                            'V': None,       # y direction mean velocity
                            'W': None,       # z direction mean velocity
                            'P': None,       # in-plane velocities (U and V components, not W)
                            'M': None,       # velocity magnitude (all three components)

                            'u': None,       # mean absolute fluctuation in U
                            'v': None,       # mean absolute fluctuation in V
                            'w': None,       # mean absolute fluctuation in W

                            'uu': None,      # turbulent energy in u (u' * u') bar
                            'vv': None,      # turbulent energy in v (v' * v') bar
                            'ww': None,      # turbulent energy in w (w' * w') bar
                            'ctke': None,    # cartesian turbulent kinetic energy (TKE)

                            'uv': None,      # reynolds stress in u/v
                            'uw': None,      # reynolds stress in u/w
                            'vw': None}      # reynolds stress in v/w

        # Build up the data set
        if v3d_paths is not None:
            self.ingest_paths(v3d_paths, min_points)  # creates constituent_vel_matrix_list from input filepath list


    def __getitem__(self, key):
        """ allows components of the instance to be accessed more simply through instance[key] """
        if key in self.vel_matrix.keys():
            return self.vel_matrix[key]

        # allows 'uv' to be returned for input key 'vu', which is nice.
        elif key[::-1] in self.vel_matrix.keys():
            return self.vel_matrix[key[::-1]]

        elif key in self.meshgrid.keys():
            return self.meshgrid[key]


    def __setitem__(self, key, value):
        """ allows components of the instance to be set more briefly """
        if key in self.vel_matrix.keys():
            self.vel_matrix[key] = value
        elif key[::-1] in self.vel_matrix.keys():
            self.vel_matrix[key[::-1]] = value
        else:
            raise AttributeError("instance does not accept __setitem__ for '{0}'".format(key))


    def ingest_paths(self, filepath_list, min_points=None):
        """
        Creates VecField3d objects from each filepath in the list
        :param filepath_list:   list of filepaths to ingest
        """

        t = Timer()

        # inherit dimensional attributes from first cartesian field and length of file list
        first_vf = VecFieldCartesian(filepath_list[0], velocity_fs=self.velocity_fs)
        self.x_set = first_vf.x_set
        self.y_set = first_vf.y_set
        self.dims = first_vf.dims + tuple([len(filepath_list)])
        self.meshgrid = first_vf.meshgrid

        # fill the dynamic data matrix with zeros
        for key in self.dynamic_set.keys():
            self.dynamic_set[key] = np.ma.zeros(self.dims)

        # populate some attributes of the dynamic set then delete the vector field instance
        for i, filepath in enumerate(filepath_list):
            cvm = VecFieldCartesian(filepath, velocity_fs=self.velocity_fs)
            self.dynamic_set['U'][:, :, i] = cvm['U']
            self.dynamic_set['V'][:, :, i] = cvm['V']
            self.dynamic_set['W'][:, :, i] = cvm['W']
            self.dynamic_set['M'][:, :, i] = (cvm['U'] ** 2 + cvm['V'] ** 2 + cvm['W'] ** 2) ** 0.5   # magnitudes
            self.dynamic_set['P'][:, :, i] = (cvm['U'] ** 2 + cvm['V'] ** 2) ** 0.5                   # in plane mags
            del cvm

        # now extract the average and fluctuating components from the dynamic data set.
        self._get_average_and_fluctuating(min_points)
        t.finish()


    def _get_average_and_fluctuating(self, min_points):
        """  populates all the components in the matrices dynamic_set and vel_matrix """

        if min_points is None:
            min_points = DEFAULT_MIN_POINTS

        # build a new mask based on whether or not more data than min_points is available (minimum point mask)
        mpm = self.dynamic_set['U'].count(axis=2) <= min_points
        self['num'] = np.ma.masked_array(self.dynamic_set['U'].count(axis=2), mask=mpm)

        # take averages
        for component in ['U', 'V', 'W', 'M', 'P']:
            self.vel_matrix[component] = masked_rms(self.dynamic_set[component], axis=2, mask=mpm)

        # find dynamic set fluctuations by subtracting out averages
        for i in range(0, self.dims[-1]):       # cant figure out fully vectorized element wise subtraction
            self.dynamic_set['u'][:, :, i] = self.dynamic_set['U'][:, :, i] - self.vel_matrix['U']
            self.dynamic_set['v'][:, :, i] = self.dynamic_set['V'][:, :, i] - self.vel_matrix['V']
            self.dynamic_set['w'][:, :, i] = self.dynamic_set['W'][:, :, i] - self.vel_matrix['W']

        # find dynamic reynolds stresses and turbulences
        for component in ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']:
            self.dynamic_set[component] = self.dynamic_set[component[0]] * self.dynamic_set[component[1]]

        # get total cartesian turbulent kinetic energy
        self.dynamic_set['ctke'] = (self.dynamic_set['uu'] ** 2 +
                                    self.dynamic_set['vv'] ** 2 +
                                    self.dynamic_set['ww'] ** 2) / 2

        # and now take the time averaged versions of each member of the dynamic set
        for component in ['u', 'v', 'w', 'uu', 'vv', 'ww', 'ctke', 'uv', 'uw', 'vw']:
            self.vel_matrix[component] = masked_rms(self.dynamic_set[component], axis=2, mask=mpm)


if __name__ == "__main__":

    directory = "../../data_test"
    paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".v3d")]

    mvf = MeanVecFieldCartesian("Station_1", paths[0:5], min_points=1)
    print mvf['U'].count()

