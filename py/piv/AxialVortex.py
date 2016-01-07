__author__ = 'Jwely'

# general imports
import cPickle
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import spectrum

# local imports
from py.piv.MeanVecFieldCartesian import MeanVecFieldCartesian
from py.utils import cart2cyl_vector, masked_rms, masked_mean, shorthand_to_tex
from py.config import *


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
        self.name_tag = name_tag

        # vortex cylindrical specific attributes
        self.core_location = (None, None)       # position of core
        self.core_index = (None, None)          # fractional index position of core
        self.velocity_fs = velocity_fs          # free stream velocity (experimental input)
        self.core_radius = None                 # distance between core and Tmax location (mm)
        self.Tmax = None                        # maximum tangential velocity
        self.Wcore = None                       # axial velocity at the core
        self.circulation_strength = None        # -flag

        # update the coordinate meshgrid
        self.meshgrid.update({"r_mesh": None,       # radial meshgrid
                              "t_mesh": None,       # tangential meshgrid (-pi to pi about right horizontal)
                              "hv_mesh": None,      # tangential meshgrid (0 near horizontal, pi/2 for vertical)
                              "t_meshd": None,      # t_mesh converted to degrees (-180 to 180)
                              "hv_meshd": None})    # hv_mesh converted to degrees (0 to 90)

        # add to the velocity matrix and flattened version
        self.mean_set.update({'R': None,  # mean radial velocity around vortex core
                              'T': None,  # mean tangential velocity around vortex core

                              'r': None,  # fluctuation in R
                              't': None,  # fluctuation in T

                              'rr': None,  # turbulent energy in r (r' * r') bar
                              'tt': None,  # turbulent energy in t (t' * t') bar

                              'rt': None,  # reynolds stress in r/t
                              'rw': None,  # reynolds stress in r/w
                              'tw': None})  # reynolds stress in t/w

        self.dynamic_set.update({'R': None,  # mean radial velocity around vortex core
                                 'T': None,  # mean tangential velocity around vortex core

                                 'r': None,  # fluctuation in R
                                 't': None,  # fluctuation in T

                                 'rr': None,  # turbulent energy in r (r' * r') bar
                                 'tt': None,  # turbulent energy in t (t' * t') bar

                                 'rt': None,  # reynolds stress in r/t
                                 'rw': None,  # reynolds stress in r/w
                                 'tw': None})  # reynolds stress in t/w


    def to_pickle(self, pickle_path, include_dynamic=False):
        """ dumps the contents of this object to a pickle """

        # delete the constituent objects with hundreds of additional matrices to reduce pkl size
        if not include_dynamic:
            del self.dynamic_set
        else:
            for key in self.dynamic_set.keys():
                if key not in DYNAMIC_INCLUDES:
                    self.dynamic_set[key] = None

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
        new_instance.characterize()
        return new_instance


    def _rrange_parser(self, r_range):
        """ parses r_range to allow arguments of core radii """

        # checks for string inputs for r_range and converts to units of mm one element at a time
        r_range_list = list(r_range)
        for i, r_element in enumerate(r_range_list):
            if isinstance(r_element, str):
                if "r" in r_element:
                    r_range_list[i] = float(r_element.replace('r', '')) * self.core_radius
                else:
                    raise Exception("Str inputs only accepted with a trailing 'r'")
        r_range = tuple(r_range_list)
        return r_range


    def _getitem_by_rt(self, component, r_range=None, t_range=None, symmetric=None):
        """
        Subsets a specific component by radius and angle theta. Defaults to subset to 50mm from
        the vortex core around the entire 360 degree range. Note that the t_range is counter clockwise
        from the first tuple entry to the second: so (-90, 90) covers the right half of the vortex, while
        (90, -90) will cover the left half of the vortex. If symmetric is set to True, the hv conversion
        of theta will be used, which ranges from 0 to 90, where 0 is near horizontal (either way) and 90
        is near vertical (either way). An input of (30, 60) will result in a set of four wedges coming out
        of each quadrant of the vortex centered about the 45 degree line, like a beach ball.

        :param component:   component to subset, all valid __getitem__ inputs will work here
        :param r_range:     tuple of (min, max) radius range in mm from the core. Supports string arguments
                            followed with a "r" to indicate units of core radii; so ('0.9r', '1.1r') will
                            set the r_range automatically to units of mm based on the calculated core radius.
        :param t_range:     tuple of (min, max) theta range in degrees about the core.
        :param symmetric:   if True, the hv_mesh is used instead of t_mesh.

        :return: a copy of the input component matrix, but with a new mask covering r,t input conditions
        """

        # set default values when None is passed
        if t_range is None:
            t_range = (-180, 180)
        if symmetric is None:
            symmetric = False

        if r_range is None:
            r_range = (0, 50)
        else:
            r_range = self._rrange_parser(r_range)

        # apply the distance mask based on the radius range
        distance_mask = np.logical_or(r_range[0] > self['r_mesh'], self['r_mesh'] > r_range[1])

        # apply the angular mask based on the theta range and symmetric criteria
        if symmetric:
            angle_mask = np.logical_or(self['hv_meshd'] < t_range[0], t_range[1] < self['hv_meshd'])

        else:
            if t_range[1] > t_range[0]:
                angle_mask = np.logical_or(self['t_meshd'] < t_range[0], t_range[1] < self['t_meshd'])
            else:
                angle_mask = np.logical_and(t_range[0] < self['t_meshd'], self['t_meshd'] < t_range[1])

        combined_mask = np.ma.mask_or(distance_mask, angle_mask)

        # take the subsets, allow a string component or a custom input numpy array. (for dynamic data)
        if isinstance(component, str):
            rt_subset_component = np.ma.masked_array(self[component], mask=combined_mask)
        elif 'numpy' in str(type(component)):
            rt_subset_component = np.ma.masked_array(component, mask=combined_mask)
        else:
            raise Exception("Cannot understand input 'compnonet' of type {0}".format(type(component)))

        return rt_subset_component


    def characterize(self, verbose=True):
        """
        Characterizes this vortex with a variety of scalar metrics such as the radius,
        the maximum tangential velocity, the maximum axial velocities, axial velocities at
        the very center of the core, etc.

        :param verbose: prints outputs
        :return characteristics_dict:
        """

        sub_t = self._getitem_by_rt('T')
        self.Tmax = np.ma.max(sub_t)
        Tmax_location = np.unravel_index(sub_t.argmax(), sub_t.shape)
        self.core_radius = self['r_mesh'][Tmax_location]
        self.Wcore = self._getitem_by_rt('W', r_range=(0, 10)).min()

        if verbose:
            message_fmt = "Core specs: radius={r:2.2f}mm, Tmax={t:2.2f}, Wmin={w:2.2f}, Vfree={vf:2.2f}"
            print(message_fmt.format(r=self.core_radius, t=self.Tmax, w=self.Wcore, vf=self.velocity_fs))

        char_dict = {"Tmax": self.Tmax,
                     "CoreRadius": self.core_radius,
                     "Wcore": self.Wcore,
                     "Vfree": self.velocity_fs}
        return char_dict


    def _find_core(self, crange=7):
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


    def _get_cylindrical_meshgrids(self, core_location_tuple=None):
        """
        Creates cylindrical meshgrids from a core location and the existing x,y meshgrids.
        These meshgrids are simply stored as attributes of the axial vortex.

        :param core_location_tuple: location of the core in (x,y) in units of millimeters
        """

        if core_location_tuple is None:
            xc, yc = self._find_core()
        else:
            xc, yc = core_location_tuple

        self.meshgrid['r_mesh'] = ((self['x_mesh'] - xc) ** 2 + (self['y_mesh'] - yc) ** 2) ** 0.5
        self.meshgrid['t_mesh'] = np.arctan2((self['y_mesh'] - yc), (self['x_mesh'] - xc))

        self.meshgrid['hv_mesh'] = abs(self.meshgrid['t_mesh'])
        left_half = self.meshgrid['hv_mesh'] > (math.pi / 2)
        self.meshgrid['hv_mesh'][left_half] = math.pi - self.meshgrid['hv_mesh'][left_half]
        self.meshgrid['t_meshd'] = self.meshgrid['t_mesh'] * 180 / math.pi      # degrees version
        self.meshgrid['hv_meshd'] = self.meshgrid['hv_mesh'] * 180 / math.pi    # degrees version


    def build_cylindrical(self, core_location_tuple=None):
        """
        Converts cartesian coordinate attributes into cylindrical attributes, and

        :param core_location_tuple:   tuple (X mm, Y mm) of actual core location on the meshgrid
        """

        self._get_cylindrical_meshgrids(core_location_tuple)

        # use the num attribute to get the minimum point mask
        mpm = self['num'].mask

        # set up empty matrices
        new_keys = ['R', 'T', 'r', 't', 'rr', 'tt', 'rt', 'rw', 'tw']
        for key in new_keys:
            self.dynamic_set[key] = np.ma.zeros(self.dims)
            self.mean_set[key] = np.ma.zeros(self.dims)

        # populate the dynamic set first with the cylindrical conversions of
        for i in range(self.dims[-1]):
            self.dynamic_set['R'][:, :, i], self.dynamic_set['T'][:, :, i] = \
                cart2cyl_vector(self.dynamic_set['U'][:, :, i],
                                self.dynamic_set['V'][:, :, i],
                                self.meshgrid['t_mesh'])

        # now average the radial and tangential components
        self.mean_set['R'] = masked_mean(self.dynamic_set['R'], axis=2, mask=mpm)
        self.mean_set['T'] = masked_mean(self.dynamic_set['T'], axis=2, mask=mpm)

        # find dynamic set fluctuations by subtracting out averages
        for i in range(0, self.dims[-1]):       # cant figure out fully vectorized element wise subtraction
            self.dynamic_set['r'][:, :, i] = self.dynamic_set['R'][:, :, i] - self.mean_set['R']
            self.dynamic_set['t'][:, :, i] = self.dynamic_set['T'][:, :, i] - self.mean_set['T']

        self.mean_set['r'] = masked_rms(self.dynamic_set['r'], axis=2, mask=mpm)
        self.mean_set['t'] = masked_rms(self.dynamic_set['t'], axis=2, mask=mpm)

        # find dynamic reynolds stresses and turbulence
        for component in ['rr', 'tt', 'rt', 'rw', 'tw']:
            self.dynamic_set[component] = self.dynamic_set[component[0]] * self.dynamic_set[component[1]]
            self.mean_set[component] = masked_rms(self.dynamic_set[component], axis=2, mask=mpm)

        # now characterize the vortex with some important but simple statistics
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


    def _get_dynamic_subsets(self, component_y, r_range=None, t_range=None, symmetric=None):
        """
        Private method for building a subset of a dynamic data set for the purposes of scatter plotting in
        the time domain. Invoked by `scatter_plot_dynamic`. Invokes _getitem_by_rt().

        :param component_y:     component to plot on the y axis
        :param r_range:         radial dimensions to subset by
        :param t_range:         angular dimensions to subset by (theta)
        :param symmetric:       bool. Subset symmetrically? see _getitem_by_rt()
        :return: Two one dimensional vectors with y,t component set in that order.
        """

        # take care of default values
        time_step = 1 / float(SAMPLING_RATE)  # convert sampling rate to time interval between points t = 1/f

        # build up a y_set, c_set, t_set, all one dimensional arrays of the same
        # length with scatter plot data and colors. allows subseting by radial and tangential dimension.
        kwargs = {"r_range": r_range, "t_range": t_range, "symmetric": symmetric}
        t_set = []
        y_sets = {"mean": [],
                  "median": [],
                  "p95": [],
                  "p05": []}

        # get the rest of the layers (component y and the time axis)
        for i in range(0, self.dims[-1]):
            ysi = self._getitem_by_rt(self.dynamic_set[component_y][:, :, i], **kwargs).flatten()
            ysi = ysi[np.invert(ysi.mask)]
            y_sets['mean'].append(np.ma.mean(ysi))
            y_sets['median'].append(np.ma.median(ysi))
            y_sets['p95'].append(np.nanpercentile(ysi, 95))
            y_sets['p05'].append(np.nanpercentile(ysi, 05))

            t_set.append(time_step * i)

        return y_sets, t_set


    def _get_vrange(self, component, low_percentile=None, high_percentile=None):
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

        if low_percentile is None:
            low_percentile = VRANGE_DEFAULT[0]
        if high_percentile is None:
            high_percentile = VRANGE_DEFAULT[1]

        comp = self._getitem_by_rt(component).astype('float')
        vmin = np.nanpercentile(comp.filled(np.nan), low_percentile)
        vmax = np.nanpercentile(comp.filled(np.nan), high_percentile)
        return vmin, vmax


    def dynamic_plot(self, component_y, title=None, y_label=None, cmap=None,
                             y_range=None, r_range=None, t_range=None, symmetric=None,
                             tight=False, figsize=None, outpath=None):
        """
        A scatter plot of component y where the x axis is locked as time.
        :param component_y:
        :param title:
        :return:
        """

        r_range = self._rrange_parser(r_range)

        print("generating dynamic plot...")
        if title is None:
            title = "{0} over time".format(shorthand_to_tex(component_y))
        if y_label is None:
            y_label = component_y
        if figsize is None:
            figsize = (14, 5)
        if cmap is None:
            cmap = SCATTER_DEFAULT_CMAP
        x_label = "Time (seconds)"

        subset_kwargs = {'r_range': r_range, 't_range': t_range, 'symmetric': symmetric}
        y_sets, t_set = self._get_dynamic_subsets(component_y, **subset_kwargs)

        # now make the figure
        fig = plt.figure(figsize=figsize, dpi=200, facecolor='w')

        # first plot
        gs = plt.GridSpec(100, 100, bottom=0.15, left=0.02, right=0.98)
        ax1 = fig.add_subplot(gs[:, 5:24])
        xlims, ylims = self._get_plot_lims(40, 66)
        plt.ylim(ylims)
        plt.xlim(xlims)

        areaplot = ax1.contourf(self['x_mesh'], self['y_mesh'],
                                self._getitem_by_rt(component_y, **subset_kwargs),
                                CONTOUR_DEFAULT_LEVELS,
                                cmap=CONTOUR_DEFAULT_CMAP)

        circ = plt.Circle(self.core_location, radius=self.core_radius, edgecolor='k',
                          linestyle=':', facecolor='none', label="Core Boundary")
        ax1.add_patch(circ)
        plt.title("Sample Area Average")
        plt.legend(loc=4)
        plt.xlabel('X position (mm)')
        plt.ylabel('Y position (mm)')

        # second dynamic plot
        ax2 = fig.add_subplot(gs[:, 30:77])
        ax2.plot(t_set, y_sets['mean'], 'k-', label='Mean')
        ax2.plot(t_set, y_sets['p05'], 'k:', label='90% of Values')
        ax2.plot(t_set, y_sets['p95'], 'k:')
        plt.grid()
        plt.legend(loc=1)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # third plot psd
        ax3 = fig.add_subplot(gs[:, 84:100])
        ax3.psd(y_sets['mean'], NFFT=64, Fs=SAMPLING_RATE)
        plt.title("log PSD")


        if outpath:
            plt.savefig(outpath)
            plt.close()
            print("saved figure to {0}".format(outpath))
        else:
            plt.show()
            plt.close()
        return


    def scatter_plot(self, component_x, component_y, component_c=None,
                     title=None, x_label=None, y_label=None, c_label=None, cmap=None,
                     x_range=None, y_range=None, r_range=None, t_range=None, symmetric=None,
                     tight=False, figsize=None, outpath=None):
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
        :param x_range:          custom x range tuple for x axis
        :param y_range:          custom y range tuple for y axis
        :param tight:           set true to squeeze the figures margins and decrease whitespace
        :param figsize:         figure size in inches at 120 dpi.
        :return:
        """

        if title is None:
            title = "{0} vs {1}".format(shorthand_to_tex(component_y), shorthand_to_tex(component_x))
        if x_label is None:
            x_label = shorthand_to_tex(component_x)
        if y_label is None:
            y_label = shorthand_to_tex(component_y)
        if c_label is None and component_c is not None:
            c_label = shorthand_to_tex(component_c)
        if figsize is None:
            figsize = (12, 6)
        if cmap is None:
            cmap = SCATTER_DEFAULT_CMAP

        x = self._getitem_by_rt(component_x, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()
        y = self._getitem_by_rt(component_y, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()

        fig = plt.figure(figsize=figsize, dpi=120, facecolor='w', edgecolor='k')
        if component_c is not None:
            c = self._getitem_by_rt(component_c, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()
            vmin, vmax = self._get_vrange(component_c)
            plt.scatter(x, y, marker='x', c=c, cmap=cmap, vmax=vmax, vmin=vmin)
            cb = plt.colorbar(orientation='vertical')
            cb.set_label(c_label)
        else:
            plt.scatter(x, y, marker='x', color='black')

        # apply manual specifications of x and y range, otherwise guess.
        if y_range is None:
            vmin, vmax = self._get_vrange(component_y, 0, 100)
            plt.ylim(vmin - 0.1, vmax * 1.1)
        else:
            plt.ylim(y_range[0], y_range[1])

        if x_range is None:
            vmin, vmax = self._get_vrange(component_x, 0, 100)
            plt.xlim(vmin - 0.1, vmax * 1.1)
        else:
            plt.xlim(x_range[0], x_range[1])

        if tight:
            plt.tight_layout()

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)

        if outpath:
            plt.savefig(outpath)
            print("saved figure to {0}".format(outpath))
            plt.close()
        else:
            plt.show()
        return


    def scatter_plot_qual(self, component_x, component_y):
        """
        Prints quick simple scatter plot of component_x vs component_y with the points colored
        according to the number of samples making up data from that point. Useful for evaluating
        trends and differentiating between real trends and potentially spurious features.

        params are exactly as scatter_plot()
        """
        self.scatter_plot(component_x, component_y, 'num', c_label="Number of Samples")


    def quiver_plot(self, title=None, outpath=None):
        """
        creates a quiver plot of the vector field
        """

        if title is None:
            title = "Quiver plot"

        plt.figure()
        plt.quiver(self['x_mesh'], self['y_mesh'],
                   self['U'], self['V'], #self['P'],
                   #cmap=cm.jet,
                   color='blue',
                   scale=400,
                   width=0.001,
                   headwidth=2,
                   headlength=1,
                   minshaft=1)

        plt.title(title)
        plt.tight_layout()
        xlims, ylims = self._get_plot_lims(50, 50)
        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.xlabel("$X$ (mm)")
        plt.ylabel("$Y$ (mm)")

        if outpath:
            plt.savefig(outpath)
            print("saved figure to {0}".format(outpath))
            plt.close()
        else:
            plt.show()


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
                       density=[len(self.x_set) / 20, len(self.y_set) / 20])

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
            print("saved figure to {0}".format(outpath))
            plt.close()
        else:
            plt.show(fig)
        return


    def contour_plot(self, component, title=None, outpath=None):
        """
        creates a contour plot of input component

        :param component:   component to plot, any member of the meshgrid or the mean_set
        :param title:       custom title to place on the contour plot
        :param outpath:     an output path to save a png file of this plot
        :return:
        """

        if title is None:
            title = self.shorthand_to_tex(component)

        fig, ax = plt.subplots()
        vmin, vmax = self._get_vrange(component)

        cf = plt.contourf(self['x_mesh'], self['y_mesh'],
                          self._getitem_by_rt(component, CONTOUR_DEFAULT_RRANGE),
                          CONTOUR_DEFAULT_LEVELS, cmap=CONTOUR_DEFAULT_CMAP, vmin=vmin, vmax=vmax)
        cf.set_clim(vmin=vmin, vmax=vmax)
        plt.colorbar(cf)
        plt.title(title)
        plt.xlabel("$X$ ($mm$)")
        plt.ylabel("$Y$ ($mm$)")

        # plot the core location for reference
        if self.core_location[0] is not None:
            ax.scatter(*self.core_location, marker='+', s=100, c='white')

        plt.tight_layout()
        xlims, ylims = self._get_plot_lims()
        plt.xlim(xlims)
        plt.ylim(ylims)

        if outpath:
            plt.savefig(outpath)
            print("saved figure to {0}".format(outpath))
            plt.close()
        else:
            plt.show()
        return


if __name__ == "__main__":


    run = 1
    directory = r"E:\Data2\Ely_May28th\Vector\{0}".format(run)
    paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".v3d")]

    small_pkl = r"C:\Users\Jeff\Desktop\Github\thesis-pivpr\pickles\Station_{0}_test_small.pkl".format(run)
    mvf = AxialVortex("Station_{0}".format(run), paths, velocity_fs=15.22)
    mvf.to_pickle(small_pkl, include_dynamic=True)

    mvf = AxialVortex().from_pickle(small_pkl)
    mvf.build_cylindrical()

    mvf.stream_plot()
    mvf.contour_plot('cte')
    mvf.contour_plot('yte')


