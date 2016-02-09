__author__ = 'Jwely'

import os
import numpy as np
import matplotlib.pyplot as plt

from py.piv.VecFieldCartesian import VecFieldCartesian
from py.utils import Timer, masked_rms, masked_mean, get_spatial_derivative
from py.config import *


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

        # attributes that will be inherited from the constituent VecField instances
        self.x_set = None
        self.y_set = None
        self.dims = (None, None, None)
        self.meshgrid = {"x_mesh": None,
                         "y_mesh": None}

        # dictionary of matrices by key symbol. Capitols are averages, lowercase are fluctuations
        self.mean_set = {'U': None,       # x direction mean velocity
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

        # the same as mean_set, except every every time slice is included for an n+1 dimensional matrix
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

        # partial spatial derivatives. the changes with respect to Z will require comparison between
        # axial vortex datasets taken at multiple stream-wise positions.
        self.derivative_set = {'dudx': None,
                               'dudy': None,
                               'dudz': None,    # special case
                               'dvdx': None,
                               'dvdy': None,
                               'dvdz': None,    # special case
                               'dwdx': None,
                               'dwdy': None,
                               'dwdz': None}    # special case

        self.equation_terms = {'turb_visc': None}  # based on simple turbulent viscosity assumption from 1877


        # Build up the data set
        if v3d_paths is not None:
            self.ingest_paths(v3d_paths, min_points)  # creates constituent_vel_matrix_list from input filepath list

    def __getitem__(self, key):
        """ allows components of the instance to be accessed more simply through instance[key] """

        if key in self.mean_set.keys():
            return self.mean_set[key]

        # allows 'uv' to be returned for input key 'vu', which is nice.
        elif key[::-1] in self.mean_set.keys():
            return self.mean_set[key[::-1]]

        elif key in self.meshgrid.keys():
            return self.meshgrid[key]

        elif key in self.derivative_set.keys():
            return self.derivative_set[key]

        elif key in self.equation_terms.keys():
            return self.equation_terms[key]

    def __setitem__(self, key, value):
        """ allows components of the instance to be set more briefly """

        if key in self.mean_set.keys():
            self.mean_set[key] = value

        elif key[::-1] in self.mean_set.keys():
            self.mean_set[key[::-1]] = value

        elif key in self.derivative_set.keys():
            self.derivative_set[key] = value

        elif key in self.equation_terms.keys():
            self.equation_terms[key] = value

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
            self.dynamic_set[key] = np.ma.zeros(self.dims, 'float')

        # populate some attributes of the dynamic set then delete the vector field instance
        for i, filepath in enumerate(filepath_list):
            cvm = VecFieldCartesian(filepath, velocity_fs=self.velocity_fs)
            self.dynamic_set['U'][:, :, i] = cvm['U']
            self.dynamic_set['V'][:, :, i] = cvm['V']
            self.dynamic_set['W'][:, :, i] = cvm['W']
            self.dynamic_set['M'][:, :, i] = (cvm['U'] ** 2 + cvm['V'] ** 2 + cvm['W'] ** 2) ** 0.5  # magnitudes
            self.dynamic_set['P'][:, :, i] = (cvm['U'] ** 2 + cvm['V'] ** 2) ** 0.5                  # in plane mags
            del cvm

        # now extract the average and fluctuating components from the dynamic data set.
        self._get_average_and_fluctuating(min_points)
        t.finish()

    def _get_average_and_fluctuating(self, min_points):
        """  populates all the components in the matrices dynamic_set and mean_set """

        if min_points is None:
            min_points = DEFAULT_MIN_POINTS

        # build a new mask based on whether or not more data than min_points is available (minimum point mask)
        mpm = self.dynamic_set['U'].count(axis=2) <= min_points
        self['num'] = np.ma.masked_array(self.dynamic_set['U'].count(axis=2), mask=mpm)

        # take averages
        for component in ['U', 'V', 'W', 'M', 'P']:
            self.mean_set[component] = masked_mean(self.dynamic_set[component], axis=2, mask=mpm)

        # find dynamic set fluctuations by subtracting out averages
        for i in range(0, self.dims[-1]):       # cant figure out fully vectorized element wise subtraction
            self.dynamic_set['u'][:, :, i] = self.dynamic_set['U'][:, :, i] - self.mean_set['U']
            self.dynamic_set['v'][:, :, i] = self.dynamic_set['V'][:, :, i] - self.mean_set['V']
            self.dynamic_set['w'][:, :, i] = self.dynamic_set['W'][:, :, i] - self.mean_set['W']

        # find dynamic reynolds stresses and turbulence
        for components in ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']:
            self.dynamic_set[components] = self.dynamic_set[components[0]] * self.dynamic_set[components[1]]

        # get total cartesian turbulent kinetic energy
        self.dynamic_set['ctke'] = 0.5 * (self.dynamic_set['uu'] + self.dynamic_set['vv'] + self.dynamic_set['ww'])

        # and now take the time averaged versions of each member of the dynamic set
        for component in ['u', 'v', 'w', 'uu', 'vv', 'ww', 'ctke', 'uv', 'uw', 'vw']:
            self.mean_set[component] = masked_rms(self.dynamic_set[component], axis=2, mask=mpm)

    def _get_spatial_derivatives(self):
        """
        This function computes 6 of the 9 spatial derivatives needed for experimental
        validation of the turbulent viscosity hypothesis. The others are left zero, which
        can be acceptable assumptions.
        """
        xmesh_m = self['x_mesh'] / 1000
        ymesh_m = self['y_mesh'] / 1000

        self.derivative_set['dudx'], self.derivative_set['dudy'] = \
            get_spatial_derivative(self['U'], xmesh_m, ymesh_m)

        self.derivative_set['dvdx'], self.derivative_set['dvdy'] = \
            get_spatial_derivative(self['V'], xmesh_m, ymesh_m)

        self.derivative_set['dwdx'], self.derivative_set['dwdy'] = \
            get_spatial_derivative(self['W'], xmesh_m, ymesh_m)

        self.derivative_set['dudz'] = self['W'] * 0
        self.derivative_set['dvdz'] = self['W'] * 0
        self.derivative_set['dwdz'] = self['W'] * 0

        return self.derivative_set

    def get_cart_turbulent_viscosity(self):
        """
        This function estimates turbulent viscosity with the terms available.
        The assumptions made here and the reasoning behind it are complex. See thesis
        document for better understanding
        :return:
        """

        # populate the spatial derivative attributes
        self._get_spatial_derivatives()

        # gather up all of our terms
        uu = self['uu']
        vv = self['vv']
        ww = self['ww']
        uw = self['uw']
        uv = self['uv']
        vw = self['vw']
        dudx = self['dudx']
        dudy = self['dudy']
        dudz = self['dudz']             # ok to assume zero
        dvdx = self['dvdx']
        dvdy = self['dvdy']
        dvdz = self['dvdz']             # ok to assume zero
        dwdx = self['dwdx']
        dwdy = self['dwdy']
        dwdz = self['dwdz'] * 0 + 1     # convert this to one so we can still see ww variation

        '''
        print 'uu', uu.max(), uu.min()
        print 'vv', vv.max(), vv.min()
        print 'ww', ww.max(), ww.min()
        print 'uv', uv.max(), uv.min()
        print 'uw', uw.max(), uw.min()
        print 'vw', vw.max(), vw.min()

        print 'dudx', dudx.max(), dudx.min()
        print 'dudy', dudy.max(), dudy.min()
        print 'dudz', dudz.max(), dudz.min()
        print 'dvdx', dvdx.max(), dvdx.min()
        print 'dvdy', dvdy.max(), dvdy.min()
        print 'dvdz', dvdz.max(), dvdz.min()
        print 'dwdx', dwdx.max(), dwdx.min()
        print 'dwdy', dwdy.max(), dwdy.min()
        print 'dwdz', dwdz.max(), dwdz.min()
        '''

        tv = - (1 / 3) * (uu / dudx + vv / dvdy + ww / dwdz) - \
             2 * ((uw / (dudz + dwdx)) + (uv / (dudy + dvdx)) + (vw / (dvdz + dwdy)))

        self['turb_visc'] = np.abs(tv)
        return tv


    def show_heatmap(self, component):
        """ prints a quick simple heads up heatmap of input component of the mean_set attribute"""
        fig, ax = plt.subplots()
        plt.pcolor(self[component])    # see __getitem__
        plt.colorbar()
        plt.title(component)
        plt.show()


if __name__ == "__main__":

    directory = "../../data_full/55"
    paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".v3d")]

    mvf = MeanVecFieldCartesian("Station_1", paths[0:20], min_points=5)
    mvf.show_heatmap('U')
    mvf.show_heatmap('u')
    mvf.show_heatmap('W')
    mvf.show_heatmap('w')
    mvf.show_heatmap('ctke')




