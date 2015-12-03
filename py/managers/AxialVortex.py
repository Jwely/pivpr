__author__ = 'Jwely'

import cPickle
import os

import matplotlib.pyplot as plt
from matplotlib import cm
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

                                'rr': None,  # turbulent energy in r (r' * r') bar
                                'tt': None,  # turbulent energy in t (t' * t') bar

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


    def find_core(self, crange=40):
        """
        Attempts to find the core near the center of the matrix. The core is found by
        searching for the minimum value of in_plane velocities within :param crange:
        index units (not mm) of the image center.
        """

        # find x and y indices of the image center
        xic = int(len(self.x_set) / 2)
        yic = int(len(self.y_set) / 2)

        # subset the in plane matrix to near the image center and find minimum there
        sub_p = self['P'][(yic - crange):(yic + crange), (xic - crange):(xic + crange)]
        sub_yi_min, sub_xi_min = np.unravel_index(sub_p.argmin(), sub_p.shape)

        # now place the location in terms of the whole image
        xi_min = xic + (sub_xi_min - crange)
        yi_min = yic + (sub_yi_min - crange)

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

        # set up empty matrices
        depth = len(self.constituent_vel_matrix_list)
        r_set = np.ma.zeros(self.dims + tuple([depth]))     # r matrix (3d)
        t_set = np.ma.zeros(self.dims + tuple([depth]))     # t matrix (3d)
        w_set = np.ma.zeros(self.dims + tuple([depth]))     # w matrix (3d)
        r_set_p = np.ma.zeros(self.dims + tuple([depth]))   # r fluctuation matrix (3d)
        t_set_p = np.ma.zeros(self.dims + tuple([depth]))   # t fluctuation matrix (3d)
        w_set_p = np.ma.zeros(self.dims + tuple([depth]))   # w fluctuation matrix (3d)

        # build the cylindrical meshgrids
        self.meshgrid['r_mesh'] = ((self['x_mesh'] - xc) ** 2 + (self['y_mesh'] - yc) ** 2) ** 0.5
        self.meshgrid['t_mesh'] = np.arctan2((self['y_mesh'] - yc), (self['x_mesh'] - xc))

        # build a 3d matrix from constituent datasets
        for i, cvm in enumerate(self.constituent_vel_matrix_list):
            r_set[:, :, i], t_set[:, :, i] = cart2cyl_vector(cvm['U'], cvm['V'], self['t_mesh'])
            w_set[:, :, i] = cvm['W']

        self['R'] = np.ma.mean(r_set, axis=2)
        self['T'] = np.ma.mean(t_set, axis=2)

        # now subtract out the averages for fluctuation measurements (time averaged)
        for i, cvm in enumerate(self.constituent_vel_matrix_list):
            r_set_p[:, :, i] = r_set[:, :, i] - self['R']
            t_set_p[:, :, i] = t_set[:, :, i] - self['T']
            w_set_p[:, :, i] = w_set[:, :, i] - self['W']

        self['r'] = np.ma.mean(abs(r_set_p), axis=2)
        self['t'] = np.ma.mean(abs(t_set_p), axis=2)

        self['rr'] = np.ma.mean(r_set_p * r_set_p, axis=2)
        self['tt'] = np.ma.mean(t_set_p * t_set_p, axis=2)

        self['rt'] = np.ma.mean(r_set_p * t_set_p, axis=2)
        self['rw'] = np.ma.mean(r_set_p * w_set_p, axis=2)
        self['tw'] = np.ma.mean(t_set_p * w_set_p, axis=2)


    def getitem_corezone(self, component, core_distance=50):
        """
        Gets a component, but subset to within :param core_distance: mm from the core.
        this is nice for finding minimums and maximums and taking statistics without including
        data at the edges of the observation are where spurious values are common.

        Note, this function returns the full component matrix, but with an updated mask

        :param component:       component to subset
        :param core_distance:   max distance (mm) to keep unmasked.
        :return:
        """

        distance_mask = self['r_mesh'] > core_distance
        core_component = np.ma.masked_array(self[component], mask=distance_mask)
        return core_component


    def _get_plot_lims(self, x_core_dist=100, y_core_dist=120):
        """
        returns the plot limits for scater plots
        :param x_core_dist:
        :param y_core_dist:
        :return:
        """
        if self.core_location[0] is None:
            self.core_location = (len(self.x_set) / 2, len(self.y_set) / 2)
            raise Warning("core location not set! using image center instead!")

        xlim = (self.core_location[0] - x_core_dist, self.core_location[0] + x_core_dist)
        ylim = (self.core_location[1] - y_core_dist, self.core_location[1] + y_core_dist)

        return xlim, ylim


    def scatter_plot2(self, component_x, component_y, title=None, x_label=None, y_label=None):
        """
        prints quick simple scatter plot of component_x vs component_y. Useful for viewing data
        as a function of distance to vortex core (R) or angle around the core (T)

        :param component_x:     component to make the X axis
        :param component_y:     component to make the Y axis
        :param component_c:     component whos value will determine the color of the dots
        :param title:
        :return:
        """

        if title is None:
            title = "{0} vs {1}".format(component_y, component_x)
        if x_label is None:
            x_label = component_x
        if y_label is None:
            y_label = component_y

        x = self[component_x].flatten()
        y = self[component_y].flatten()
        c = self['num'].flatten()

        fig, ax = plt.subplots()

        plt.scatter(x, y, marker='x', c=c, cmap=cm.bone_r)
        cb = plt.colorbar(orientation='horizontal')
        cb.set_label("Quality of point (N good samples)")

        vmin, vmax = self._get_vrange(component_y)

        plt.ylim(vmin - 0.1, vmax * 2)
        plt.tight_layout(pad=2)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.show()


    def scatter_plot(self, component_x, component_y, component_c=None, title=None,
                         x_label=None, y_label=None, c_label=None, cmap=cm.hsv):
        """
        prints quick simple scatter plot of component_x vs component_y. Useful for viewing data
        as a function of distance to vortex core (R) or angle around the core (T)

        :param component_x:     component to make the X axis
        :param component_y:     component to make the Y axis
        :param component_c:     component who's value will determine the color of the dots
        :param title:
        :return:
        """

        if title is None:
            title = "{0} vs {1}".format(component_y, component_x)
        if x_label is None:
            x_label = component_x
        if y_label is None:
            y_label = component_y

        x = self[component_x].flatten()
        y = self[component_y].flatten()

        fig, ax = plt.subplots()
        if component_c is not None:
            c = self[component_c].flatten()
            vmin, vmax = self._get_vrange(component_c)
            plt.scatter(x, y, marker='x', c=c, cmap=cmap, vmax=vmax, vmin=vmin)
            cb = plt.colorbar(orientation='horizontal')
            cb.set_label(c_label)
        else:
            plt.scatter(x, y, marker='x', color='black')

        vmin, vmax = self._get_vrange(component_y)

        plt.ylim(vmin - 0.1, vmax * 2)
        plt.tight_layout(pad=2)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.show()


    def stream_plot(self, title=None):
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
                plt.scatter(*self.core_location, marker='+', s=100, c='white')

        # show the figure
        plt.show(fig)


    def contour_plot(self, components, titles=None, shape=None):

        if titles is None:
            titles = components

        # if there is just one component to plot
        if isinstance(components, str):
            self._single_contour_plot(components, titles)

        elif isinstance(components, list):
            self._multi_contour_plot(components, titles, shape)




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

    mvf.stream_plot()
    mvf.contour_plot('cte')
    mvf.contour_plot('yte')


