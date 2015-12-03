__author__ = 'Jwely'

import cPickle
import os

import matplotlib.pyplot as plt
from matplotlib import cm
import pylab

import numpy as np

from py.managers.MeanVecFieldCartesian import MeanVecFieldCartesian
from py.utils.cart2cyl_vector import cart2cyl_vector


class AxialVortex(MeanVecFieldCartesian):

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

        :param name_tag:        unique string name tag for this data set
        :param v3d_paths:       list of filepaths to v3d files
        :param velocity_fs:     free stream velocity (meters/second)
        :return:
        """

        # invoke the parent class init
        MeanVecFieldCartesian.__init__(self, name_tag, v3d_paths, velocity_fs)
        name_tag = name_tag

        # add vortex cylindrical specific attributes
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


    def to_pickle(self, pickle_path, include_dynamic=False):
        """ dumps the contents of this object to a pickle """

        # delete the constituent objects with hundreds of additional matrices to reduce pkl size
        if not include_dynamic:
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


    def find_core(self, range=40):
        """
        Attempts to find the core near the center of the matrix. The core is found by
        searching for the minimum value of in_plane velocities within :param range:
        index units (not mm) of the image center.
        """

        # find x and y indices of the image center
        xic = int(len(self.x_set) / 2)
        yic = int(len(self.y_set) / 2)

        # subset the in plane matrix to near the image center and find minimum there
        sub_p = self['P'][(yic - range):(yic + range), (xic - range):(xic + range)]
        sub_yi_min, sub_xi_min = np.unravel_index(sub_p.argmin(), sub_p.shape)

        # now place the location in terms of the whole image
        xi_min = xic + (sub_xi_min - range)
        yi_min = yic + (sub_yi_min - range)

        # subset again, in the immediate core zone to interpolate a "true" core position
        cz = self['P'][(yi_min - 1):(yi_min + 2), (xi_min - 1):(xi_min + 2)]
        cz_x_mesh = self['x_mesh'][(yi_min - 1):(yi_min + 2), (xi_min - 1):(xi_min + 2)]
        cz_y_mesh = self['y_mesh'][(yi_min - 1):(yi_min + 2), (xi_min - 1):(xi_min + 2)]

        # just take an inverse in-plane velocity weighted average of the meshgrids
        xc = np.sum((1 / cz) * cz_x_mesh) / np.sum(1 / cz)  # x coordinate of core axis
        yc = np.sum((1 / cz) * cz_y_mesh) / np.sum(1 / cz)  # y coordinate of core axis

        self.core_location = (xc, yc)
        return self.core_location


    def build_cylindrical(self, core_location_tuple=None):
        """
        Converts cartesian coordinate attributes into cylindrical attributes, and

        :param core_location_tuple:   tuple (X mm, Y mm) of actual core location on the meshgrid
        """

        if core_location_tuple is not None:
            xc, yc = core_location_tuple
        else:
            xc, yc = self.find_core()

        # build the cylindrical meshgrid
        self.meshgrid['r_mesh'] = ((self['x_mesh'] - xc) ** 2 + (self['y_mesh'] - yc) ** 2) ** 0.5
        self.meshgrid['t_mesh'] = np.arctan2((self['y_mesh'] - yc), (self['x_mesh'] - xc))

        self['R'], self['T'] = cart2cyl_vector(self['U'], self['V'], self['t_mesh'])
        self['r'], self['t'] = cart2cyl_vector(self['u'], self['v'], self['t_mesh'])

        self['rr'] = self['r'] ** 2.0
        self['tt'] = self['t'] ** 2.0
        self['yte'] = self['rr'] + self['tt'] + self['ww']

        self['rt'] = self['r'] * self['t']
        self['rw'] = self['r'] * self['w']
        self['tw'] = self['t'] * self['w']
        self['yrs'] = self['rt'] + self['rw'] + self['tw']


    def set_plot_range(self, x_extent=None, y_extent=None):
        """
        accepts kwargs for x and y extents
        :param x_extent:
        :param y_extent:
        :return:
        """


    def get_scatter_plot(self, component_y, component_x, title=None):
        """
        prints quick simple scatter plot of component_x vs component_y. Useful for viewing data
        as a function of distance to vortex core (R) or angle around the core (T) """

        if title is None:
            title = "{0} vs {1}".format(component_x, component_y)

        fig, ax = plt.subplots()




    def get_stream_plot(self, title=None):
        """
        Renders a stream plot of the data to the screen.
        :return:
        """

        if title is None:
            title = "Stream: colored by in-plane velocities"

        fig, ax = plt.subplots()
        plt.streamplot(self['x_mesh'],
                       self['y_mesh'],
                       self['U'],
                       self['V'],
                       color=self['P'],
                       arrowstyle='->',
                       arrowsize=1,
                       density=[len(self.x_set) / 20, len(self.y_set) / 20],
                       )

        plt.colorbar()
        plt.title(title)

        # plot the core location for reference
        if self.core_location[0] is not None:
            ax.scatter(*self.core_location, marker='+', s=200, c='black')

        plt.show(fig)


    def _get_vrange(self, component, low_percentile=3, high_percentile=97):
        """
        Gets the percentile range for color bar scaling on a given component matrix. Used
        to ensure spurious high and low values do not over stretch the color ramps on plots.
        The more noise in the data the further form 0 and 100 respectively these percentiles
        must be.

        :param component:           component to get range for
        :param low_percentile:      low percentile value marking coolest color
        :param high_percentile:     high percentile value marking warmest color
        :return:
        """
        vmin = np.percentile(self[component], low_percentile)
        vmax = np.percentile(self[component], high_percentile)
        return vmin, vmax


    def _single_contour_plot(self, component, title=None):
        """ Handles the instances in which only a single contour plot is desired """
        fig, ax = plt.subplots()
        vmin, vmax = self._get_vrange(component)

        cf = plt.contourf(self['x_mesh'], self['y_mesh'], self[component], 256,
                          cmap=cm.jet, vmin=vmin, vmax=vmax)
        cf.set_clim(vmin=vmin, vmax=vmax)
        plt.colorbar(cf)
        plt.title(title)
        plt.xlabel("X position (mm)")
        plt.ylabel("Y position (mm)")

        # plot the core location for reference
        if self.core_location[0] is not None:
            ax.scatter(*self.core_location, marker='+', s=100, c='white')
        plt.show()


    def _multi_contour_plot(self, components, titles=None, shape=None):
        """
        Manages multiple user plots.
        :param components:      list of components
        :param titles:          list of titles
        :param shape:           tuple of (nrows, ncols) for figure tiling matrix
        """

        # determine the shape of the subplots
        if shape is None:
            nrows = 1
            ncols = len(components)
        else:
            nrows, ncols = shape

        # create the figure
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols)

        # add each component plot
        for i, component in enumerate(components):

            plt.subplot(nrows, ncols, i + 1)
            vmin, vmax = self._get_vrange(component)

            cf = plt.contourf(self['x_mesh'], self['y_mesh'], self[component], 64,
                              cmap=cm.jet, vmin=vmin, vmax=vmax)

            plt.colorbar(cf)
            plt.xlabel("X position (mm)")
            plt.ylabel("Y position (mm)")
            plt.title(titles[i])

            # plot the core location for reference
            if self.core_location[0] is not None:
                ax[i] = plt.scatter(*self.core_location, marker='+', s=100, c='white')

        # show the figure
        plt.show(fig)


    def get_contour_plot(self, components, titles=None):

        if titles is None:
            titles = components

        # if there is just one component to plot
        if isinstance(components, str):
            self._single_contour_plot(components, titles)

        elif isinstance(components, list):
            self._multi_contour_plot(components, titles)


if __name__ == "__main__":


    run = 1
    directory = r"E:\Data2\Ely_May28th\Vector\{0}".format(run)
    paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".v3d")]

    small_pkl = r"C:\Users\Jeff\Desktop\Github\thesis-pivpr\pickles\Station_{0}_test_small.pkl".format(run)
    mvf = AxialVortex("Station_{0}".format(run), paths, velocity_fs=15.22)
    mvf.to_pickle(small_pkl, include_dynamic=True)

    mvf = AxialVortex().from_pickle(small_pkl)
    mvf.find_core()
    mvf.build_cylindrical()

    mvf.get_stream_plot()
    mvf.get_contour_plot('cte')
    mvf.get_contour_plot('yte')


