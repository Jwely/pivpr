__author__ = 'Jwely'

import cPickle
import os

import numpy as np

from MeanVecFieldCartesian import MeanVecFieldCartesian
from py.utils.cart2cyl_vector import cart2cyl_vector


class MeanVecFieldCylindrical(MeanVecFieldCartesian):

    def __init__(self, name_tag=None, v3d_paths=None, velocity_fs=None):
        """
        Built to extend the cartesian version of this class. Since all PIV data is
        reasonably always going to be taken raw in cartesian coordinates, there is no
        native cylindrical class for mean vector fields.

        This class allows the precise definition of an axis location in (x,y) space and transforms
        all (x,y) data into (r,t) data. The z coordinate should be the streamwise coordinate.

        X = span wise, right positive?
        Y = height wise, up positive?
        Z = stream wise, downstream positive?
        R = radial outwards, always positive.
        T = tangential, clockwise positive?

        :param v3d_paths:       list of filepaths to v3d files
        :param name_tag:        unique string name tag for this data set
        :return:
        """

        # invoke the parent class init
        MeanVecFieldCartesian.__init__(self, name_tag, v3d_paths, velocity_fs)

        # add cylindrical specific attributes
        self.core_location = (None, None)       # position of core
        self.core_index = (None, None)          # fractional index position of core

        self.meshgrid.update({"r_mesh": None,   # radial meshgrid
                              "t_mesh": None})  # tangential meshgrid

        # add to the velocity matrix and flattened version
        self.vel_matrix.update({'R': None,  # mean radial velocity around vortex core
                                'T': None,  # mean tangential velocity around vortex core

                                'r': None,  # fluctuation in R
                                't': None,  # fluctuation in T

                                'rr': None,  # turbulent energy in radial
                                'tt': None,  # turbulent energy in tangential
                                'yte': None,  # total cylindrical turbulent energy

                                'rt': None,  # reynolds stress in r/t
                                'rw': None,  # reynolds stress in r/w
                                'tw': None,  # reynolds stress in t/w
                                'yrs': None})  # total cylindrical reynolds stress


    def to_pickle(self, pickle_path, reduce_memory=False):
        """ dumps the contents of this object to a pickle """

        # delete the constituent objects with hundreds of additional matrices to reduce pkl size
        if reduce_memory:
            del self.constituent_vel_matrix_list
            self.constituent_vel_matrix_list = None

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


    def find_core(self, range=20):
        """
        Attempts to find the core near the center of the matrix. The core is found by searching
        for the minimum value of in_plane velocities within :param range: mm of the image center.
        :return:
        """

        # find x and y indices of the image center
        xic = len(self.x_set) / 2
        yic = len(self.y_set) / 2

        # subset the in plane matrix to near the image center and find minimum there
        sub_p = self['P'][(xic - range):(xic + range), (yic - range):(yic + range)]
        sub_xi_min, sub_yi_min = np.unravel_index(sub_p.argmin(), sub_p.shape)

        # now place the location in terms of the whole image
        xi_min = sub_xi_min + xic
        yi_min = sub_yi_min + yic

        print self.x_set[xi_min], self.y_set[yi_min]



    def build_cylindrical(self, core_location_tuple):
        """
        Converts cartesian coordinate attributes into cylindrical attributes, and

        :param core_location_tuple:   tuple (X mm, Y mm) of actual core location on the meshgrid
        """

        xc, yc = core_location_tuple

        # build the cylindrical meshgrid
        self.meshgrid['r_mesh'] = ((self['x_mesh'] - xc) ** 2 + (self['y_mesh'] - yc) ** 2) ** 0.5
        self.meshgrid['t_mesh'] = np.arctan2((self['y_mesh'] - yc), (self['x_mesh'] - xc))

        self['R'], self['T'] = cart2cyl_vector(self['U'], self['V'], self['t_mesh'])
        self['r'], self['t'] = cart2cyl_vector(self['u'], self['v'], self['t_mesh'])

        self['rr'] = self['r'] ** 2.0
        self['tt'] = self['t'] ** 2.0
        self['yte'] = self['r'] + self['t'] + self['w']

        self['rt'] = self['r'] * self['t']
        self['rw'] = self['r'] * self['w']
        self['tw'] = self['t'] * self['w']
        self['yrs'] = self['rt'] + self['rw'] + self['tw']


    def show_scatter(self, component_y, component_x):
        """
        prints quick simple scatter plot of component_x vs component_y. Useful for viewing data
        as a function of distance to vortex core (R) or angle around the core (T) """

        # will need to 1-dimensionalize the vel_matrix data for plotting purposes.
        pass


if __name__ == "__main__":

    run = 1

    directory = r"E:\Data2\Ely_May28th\Vector\{0}".format(run)
    paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".v3d")]

    small_pkl = r"C:\Users\Jeff\Desktop\Github\thesis-pivpr\pickles\Station_{0}_test_small.pkl".format(run)
    #mvf = MeanVecFieldCylindrical("Station_{0}".format(run), paths, velocity_fs=15.22)
    #mvf.to_pickle(small_pkl, reduce_memory=True)

    mvf = MeanVecFieldCylindrical().from_pickle(small_pkl)
    mvf.find_core()
    mvf.show_contour('P')
    #mvf.build_cylindrical((73.5214, 43.4737))

    #mvf.show_stream()
    #mvf.show_contour('t_mesh')

