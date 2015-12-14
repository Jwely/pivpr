__author__ = 'Jwely'

import cPickle
import os

import matplotlib.pyplot as plt
from matplotlib import cm

# this prevents plt.tight_layout() from crowding axis labels off the edges of the plot.
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

import numpy as np

from py.manager.MeanVecFieldCartesian import MeanVecFieldCartesian
from py.utils.cart2cyl_vector import cart2cyl_vector


class AxialVortex(MeanVecFieldCartesian):

    def __init__(self, name_tag=None, v3d_paths=None, velocity_fs=None, min_points=20):
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
        MeanVecFieldCartesian.__init__(self, name_tag=name_tag, v3d_paths=v3d_paths,
                                       velocity_fs=velocity_fs, min_points=min_points)
        name_tag = name_tag

        # vortex cylindrical specific attributes
        self.core_location = (None, None)       # position of core
        self.core_index = (None, None)          # fractional index position of core
        self.velocity_fs = velocity_fs          # free stream velocity (experimental input)
        self.core_radius = None                 # distance between core and Tmax location (mm)
        self.Tmax = None                        # maximum tangential velocity
        self.Wcore = None                       # axial velocity at the core
        self.circulation_strength = None        # -flag

        # update the coordinate meshgrid
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


    def from_pickle(self, pickle_path):
        """ loads previous saved state from a .pkl file and returns a MeanVecFieldCartesian instance """

        with open(pickle_path, 'rb') as f:
            new_instance = cPickle.load(f)

        print("loaded pkl from {0}".format(pickle_path))
        new_instance.characterize()
        return new_instance


    def _getitem_corezone(self, component, core_distance=50):
        """
        Gets a component, but subset to within :param core_distance: mm from the core.
        this is nice for finding minimums and maximums and taking statistics without including
        data at the edges of the observation are where spurious values are common.

        Note, this function returns the full component matrix, but with an updated mask

        :param component:       component to subset
        :param core_distance:   max distance (mm) to keep unmasked.
        :return:
        """

        if self[component] is None:
            raise Exception("cannot find component '{0}'".format(component))

        distance_mask = self['r_mesh'] > core_distance
        core_component = np.ma.masked_array(self[component], mask=distance_mask)
        return core_component


    def characterize(self, verbose=True):
        """
        Characterizes this vortex with a variety of scalar metrics such as the radius,
        the maximum tangential velocity, the maximum axial velocities, axial velocities at
        the very center of the core, etc.

        :param verbose: prints outputs
        :return characteristics_dict:
        """

        sub_t = self._getitem_corezone('T')
        self.Tmax = np.ma.max(sub_t)
        Tmax_location = np.unravel_index(sub_t.argmax(), sub_t.shape)
        self.core_radius = self['r_mesh'][Tmax_location]
        self.Wcore = self._getitem_corezone('W', 10).min()

        if verbose:
            message_fmt = "Core specs: radius={r:2.2f}mm, Tmax={t:2.2f}, Wmin={w:2.2f}, Vfree={vf:2.2f}"
            print(message_fmt.format(r=self.core_radius, t=self.Tmax, w=self.Wcore, vf=self.velocity_fs))

        char_dict = {"Tmax": self.Tmax,
                     "CoreRadius": self.core_radius,
                     "Wcore": self.Wcore,
                     "Vfree": self.velocity_fs}
        return char_dict


    def _find_core(self, crange=20):
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
            xc, yc = self._find_core()

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

        # now characterize the vortex
        characteristics = self.characterize(verbose=True)
        return characteristics


    def _get_plot_lims(self, x_core_dist=100, y_core_dist=100):
        """
        returns the plot extents based on user defined distance to core. If the
        user input distances cause areas outside the image extent to be displayed, they
        are automatically trimmed not to exceed the image extents.

        :param x_core_dist:  maximum x distance from core to show on plot
        :param y_core_dist:  maximum y distance from core to show on plot
        :return:
        """
        if self.core_location[0] is None:
            self.core_location = (len(self.x_set) / 2, len(self.y_set) / 2)
            raise Warning("core location not set! using image center instead!")

        # use either the input distance away from the core, or the image extent, whichever is limiting
        xlim_low = max([self.core_location[0] - x_core_dist, min(self.x_set)])
        xlim_high = min([self.core_location[0] + x_core_dist, max(self.x_set)])
        ylim_low = max([self.core_location[1] - y_core_dist, min(self.y_set)])
        ylim_high = min([self.core_location[1] + y_core_dist, max(self.y_set)])

        xlim = (xlim_low, xlim_high)
        ylim = (ylim_low, ylim_high)

        return xlim, ylim


    def scatter_plot_qual(self, component_x, component_y, title=None, x_label=None, y_label=None):
        """
        prints quick simple scatter plot of component_x vs component_y with the points colored
        according to the number of samples making up data from that point. Useful for evaluating
        trends and differentiating between real trends and potentially spurious features.

        :param component_x:     component to make the X axis
        :param component_y:     component to make the Y axis
        :param title:           custom title
        :param x_label:         custom x axis label
        :param y_label:         custom y axis label
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
        plt.tight_layout()
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.show()


    def scatter_plot(self, component_x, component_y, component_c=None, title=None,
                         x_label=None, y_label=None, c_label=None, cmap=cm.hsv,
                         xrange=None, yrange=None, tight=False, figsize=None, outpath=None):
        """
        prints quick simple scatter plot of component_x vs component_y. Useful for viewing data
        as a function of distance to vortex core (R) or angle around the core (T)

        :param component_x:     component to make the X axis
        :param component_y:     component to make the Y axis
        :param component_c:     component who's value will determine the color of the dots
        :param title:           custom title for plot
        :param x_label:         custom x axis label
        :param y_label:         custom y axis label
        :param c_label:         custom color bar label
        :param cmap:            custom colormap for color bar
        :param xrange:          custom x range tuple for x axis
        :param yrange:          custom y range tuple for y axis
        :param tight:           set true to squeeze the figures margins and decrease whitespace
        :param figsize:         figure size in inches at 120 dpi.
        :return:
        """

        if title is None:
            title = "{0} vs {1}".format(component_y, component_x)
        if x_label is None:
            x_label = component_x
        if y_label is None:
            y_label = component_y
        if figsize is None:
            figsize = (12, 6)

        x = self[component_x].flatten()
        y = self[component_y].flatten()


        fig = plt.figure(figsize=figsize, dpi=120, facecolor='w', edgecolor='k')
        if component_c is not None:
            c = self[component_c].flatten()
            vmin, vmax = self._get_vrange(component_c)
            plt.scatter(x, y, marker='x', c=c, cmap=cmap, vmax=vmax, vmin=vmin)
            cb = plt.colorbar(orientation='vertical')
            cb.set_label(c_label)
        else:
            plt.scatter(x, y, marker='x', color='black')


        # apply manual specifications of x and y range, otherwise guess.
        if yrange is None:
            vmin, vmax = self._get_vrange(component_y, 0, 100)
            plt.ylim(vmin - 0.1, vmax * 1.1)
        else:
            plt.ylim(yrange[0], yrange[1])

        if xrange is None:
            vmin, vmax = self._get_vrange(component_x, 0, 100)
            plt.xlim(vmin - 0.1, vmax * 1.1)
        else:
            plt.xlim(xrange[0], xrange[1])

        if tight:
            plt.tight_layout()

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)

        if outpath:
            plt.savefig(outpath)
        else:
            plt.show()
        return



    def stream_plot(self, title=None, outpath=None):
        """
        Renders a stream plot of the data to the screen.
        :param title:   A custom title for the stream plot
        """

        if title is None:
            title = "Stream: colored by in-plane velocities"

        fig, ax = plt.subplots()
        plt.streamplot(self['x_mesh'], self['y_mesh'], self['U'], self['V'],
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

        plt.tight_layout()
        xlims, ylims = self._get_plot_lims()
        plt.xlim(xlims)
        plt.ylim(ylims)
        if outpath:
            plt.savefig(outpath)
        else:
            plt.show(fig)
        return




    def _get_vrange(self, component, low_percentile=1, high_percentile=99.9):
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
        comp = self._getitem_corezone(component).astype('float')
        vmin = np.nanpercentile(comp.filled(np.nan), low_percentile)
        vmax = np.nanpercentile(comp.filled(np.nan), high_percentile)
        return vmin, vmax


    def _single_contour_plot(self, component, title=None, outpath=None):
        """ Handles the instances in which only a single contour plot is desired """
        fig, ax = plt.subplots()
        vmin, vmax = self._get_vrange(component)

        cf = plt.contourf(self['x_mesh'], self['y_mesh'],
                          self._getitem_corezone(component, core_distance=100),
                          512,
                          cmap=cm.jet, vmin=vmin, vmax=vmax)
        cf.set_clim(vmin=vmin, vmax=vmax)
        plt.colorbar(cf)
        plt.title(title)
        plt.xlabel("X position (mm)")
        plt.ylabel("Y position (mm)")

        # plot the core location for reference
        if self.core_location[0] is not None:
            ax.scatter(*self.core_location, marker='+', s=100, c='white')

        plt.tight_layout()
        xlims, ylims = self._get_plot_lims()
        plt.xlim(xlims)
        plt.ylim(ylims)

        if outpath:
            plt.savefig(outpath)
        else:
            plt.show()
        return


    def _multi_contour_plot(self, components, titles=None, shape=None, outpath=None):
        """ Manages multiple user plots """

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

            cf = plt.contourf(self['x_mesh'], self['y_mesh'], self[component], 256,
                              cmap=cm.jet, vmin=vmin, vmax=vmax)

            plt.colorbar(cf)
            plt.xlabel("X position (mm)")
            plt.ylabel("Y position (mm)")
            plt.title(titles[i])

            # plot the core location for reference
            if self.core_location[0] is not None:
                plt.scatter(*self.core_location, marker='+', s=100, c='white')

            xlims, ylims = self._get_plot_lims()
            plt.xlim(xlims)
            plt.ylim(ylims)

        # show the figure
        plt.subplots_adjust(left=0.03, right=0.97, wspace=0.1, hspace=0.2)

        if outpath:
            plt.savefig(outpath)
        else:
            plt.show(fig)
        return


    def contour_plot(self, components, titles=None, shape=None):
        """
        Creates a contour plot, accepts multiple plots in a subfigure according to shape layout

        :param components:  string or list of all components to plot
        :param titles:      string or list of custom titles to plot
        :param shape:       custom shape layout, defaults to multiple plots horizontally
        """

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
    mvf._find_core()
    mvf.build_cylindrical()

    mvf.stream_plot()
    mvf.contour_plot('cte')
    mvf.contour_plot('yte')


