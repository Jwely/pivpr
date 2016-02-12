__author__ = 'Jwely'

# general imports
import cPickle
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# package imports
from shorthand_to_tex import shorthand_to_tex
from py.piv.MeanVecFieldCartesian import MeanVecFieldCartesian
from py.vortex_theory import AshVortex, LambOseenVortex, RankineVortex, BurnhamHallockVortex
from py.utils import cart2cyl_vector, masked_rms, masked_mean, smooth_filt, get_spatial_derivative, dbz
from py.config import *


class AxialVortex(MeanVecFieldCartesian):
    def __init__(self, name_tag=None, v3d_paths=None, velocity_fs=None, z_location=None, min_points=20):
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
        self.z_location = z_location        # the position of this vortex downstream. in mm
        self.core_location = (None, None)   # position of core
        self.core_index = (None, None)      # fractional index position of core
        self.velocity_fs = velocity_fs      # free stream velocity (experimental input)
        self.core_radius = None             # distance between core and Tmax location (mm)
        self.Tmax = None                    # maximum tangential velocity
        self.Wcore = None                   # axial velocity at the core
        self.Wmean = None                   # average axial velocity within the good data range
        self.circulation_strength = None    # -flag

        self.char_dict = None               # dictionary for characteristics (used in outputing)

        # update the coordinate meshgrid
        self.meshgrid.update({"r_mesh": None,  # radial meshgrid
                              "t_mesh": None,  # tangential meshgrid (-pi to pi about right horizontal)
                              "hv_mesh": None,  # tangential meshgrid (0 near horizontal, pi/2 for vertical)
                              "t_meshd": None,  # t_mesh converted to degrees (-180 to 180)
                              "hv_meshd": None,  # hv_mesh converted to degrees (0 to 90)
                              })


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

        # standard cylindrical derivatives which aren't already covered in parent class
        self.derivative_set.update({'drdr': None,
                                    'drdt': None,
                                    'drdz': None,    # special case
                                    'dtdr': None,
                                    'dtdt': None,
                                    'dtdz': None,    # special case
                                    'dwdr': None,
                                    'dwdt': None})

        self.derivative_set.update({'sr_cx1': None,   # first cylindrical strain rate term approx on X axis
                                    'sr_cx2': None,   # second cylindrical strain rate term approx on X axis
                                    'sr_cx3': None,   # third cylindrical strain rate term approx on X axis
                                    'sr_cx4': None,   # fourth cylindrical strain rate term approx on X axis
                                    'sr_cx_tot': None,  # the sum of all cylindrical strain rate terms approx on X axis
                                    })

        self.equation_terms.update({'turb_visc_reynolds': None,  # for turbulent viscosity by pressure relaxation calc
                                    'turb_visc_vel_grad': None,  # for turbulent viscosity by pressure relaxation calc
                                    'turb_visc_ettap': None,     # the pressure relaxation term
                                    'turb_visc_total': None,     # for turbulent viscosity by pressure relaxation calc
                                    'turb_visc_ratio': None,     # the ratio between total and classical turb_visc
                                    'dPdr': None,                # The calculated radial pressure gradient
                                    })


    def to_pickle(self, pickle_path, include_dynamic=False):
        """ dumps the contents of this object to a pickle """

        # delete the constituent objects with hundreds of additional matrices to reduce pkl size
        temp_dynamic = self.dynamic_set    # save a copy of the unmodified dynamic set
        if not include_dynamic:
            del self.dynamic_set
        else:
            for key in self.dynamic_set.keys():
                if key not in DYNAMIC_INCLUDES:
                    self.dynamic_set[key] = None

        # create the directory and write the pkl file.
        pkl_dir = os.path.dirname(pickle_path)
        if pkl_dir != "" and not os.path.exists(pkl_dir):
            os.mkdir(os.path.dirname(pickle_path))

        with open(pickle_path, 'wb+') as f:
            cPickle.dump(self, f)
            print("Saved to {0}".format(pickle_path))
            
        # restore the dynamic data to this instance in case it needs to be used subsequently
            self.dynamic_set = temp_dynamic
            del temp_dynamic

    @staticmethod
    def from_pickle(pickle_path):
        """ loads previous saved state from a .pkl file and returns a MeanVecFieldCartesian instance """

        with open(pickle_path, 'rb') as f:
            new_instance = cPickle.load(f)

        print("loaded pkl from {0}".format(pickle_path))
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
            symmetric = False

        if r_range is None:
            r_range = (0, 70)
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

        # take the subsets, allow a string component or a custom input np array. (for dynamic data)
        if isinstance(component, str):
            if self[component] is None:
                raise Exception("Component {0} is None!".format(component))
            else:
                rt_subset_component = np.ma.masked_array(self[component], mask=combined_mask)
        elif 'numpy' in str(type(component)):
            rt_subset_component = np.ma.masked_array(component, mask=combined_mask)
        else:
            raise Exception("Cannot understand input component of type {0}".format(type(component)))

        return rt_subset_component

    def _get_rcore_tmax(self):
        # get the smoothed curvesmooth the curve
        r_avg, t_avg = self.get_smoothed_line('r_mesh', 'T', 1e4, 100, 50, 50)

        # chop of back 5% to remove the smoothing distortion
        r_avg = r_avg[0:-int(len(r_avg) / 20)]
        t_avg = t_avg[0:-int(len(t_avg) / 20)]
        t_max = t_avg.max()
        r_core = r_avg[np.unravel_index(np.ma.argmax(t_avg), t_avg.shape)]
        return r_core, t_max

    def _guess_core(self, crange=20):
        """
        creates an initial guess at the location of the vortex core based on in plane velocities
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
        core_location = (xc, yc)
        return core_location

    def find_core(self, crange=20, core_radius_range=(10.0, 25.0), vtheta_max_range=(3.0, 10.0), verbose=True):
        """
        Attempts to find the core near the center of the matrix by restricting the search area
        iteratively. characterizes the vortex once a core is found.

        :param crange:              the number of grid points away from geometric center to search
                                    for a minimum in plane velocity.
        :param core_radius_range:   the range in mm to look for a vtehta_max, a core radius outside this
                                    range will be considered erroneous.
        :param vtheta_max_range:    the range of acceptable vtheta_max values. values outside this
                                    # range will be considered erroneous
        """
        core_radius_good = False
        t_max_good = False
        max_tries = 10
        tries = 0

        while not core_radius_good and not t_max_good or tries > max_tries:
            core_location = self._guess_core(crange=crange)
            self._get_cylindrical_meshgrids(core_location)
            self._build_cylindrical(core_location)
            r_core, t_max = self._get_rcore_tmax()

            if core_radius_range[0] < r_core < core_radius_range[1]:
                core_radius_good = True
            if vtheta_max_range[0] < t_max < vtheta_max_range[1]:
                t_max_good = True

            crange *= 0.8
            tries += 1

        print("found core with crange={0}".format(crange))
        self.Tmax = t_max
        self.core_radius = r_core
        self.core_location = core_location
        self.Wcore = self._getitem_by_rt('W', r_range=(0, 10)).min()
        self.Wmean = self._getitem_by_rt('W', r_range=(0, 80)).mean()

        # if free stream velocity isn't specified, infer it from mean axial velocity.
        if self.velocity_fs is None:
            self.velocity_fs = self.Wmean

        if verbose:
            message_fmt = "Core specs: r={r:2.2f}mm, Tmax={t:2.2f}, Wmin={w:2.2f}, Vfree={vf:2.2f}"
            print(message_fmt.format(r=self.core_radius, t=self.Tmax, w=self.Wcore, vf=self.velocity_fs))

        char_dict = {"T_max": self.Tmax,
                     "r_mesh_core": self.core_radius,
                     "W_core": self.Wcore,
                     "W_mean": self.Wmean,
                     "velocity_free_stream": self.velocity_fs,
                     "core_location": self.core_location}
        self.char_dict = char_dict
        return char_dict

    def _get_cylindrical_meshgrids(self, core_location_tuple):
        """
        Creates cylindrical meshgrids from a core location and the existing x,y meshgrids.
        These meshgrids are simply stored as attributes of the axial vortex.

        :param core_location_tuple: location of the core in (x,y) in units of millimeters
        """

        xc, yc = core_location_tuple

        self.meshgrid['r_mesh'] = ((self['x_mesh'] - xc) ** 2 + (self['y_mesh'] - yc) ** 2) ** 0.5
        self.meshgrid['t_mesh'] = np.arctan2((self['y_mesh'] - yc), (self['x_mesh'] - xc))

        self.meshgrid['hv_mesh'] = abs(self.meshgrid['t_mesh'])
        left_half = self.meshgrid['hv_mesh'] > (math.pi / 2)
        self.meshgrid['hv_mesh'][left_half] = math.pi - self.meshgrid['hv_mesh'][left_half]
        self.meshgrid['t_meshd'] = self.meshgrid['t_mesh'] * 180 / math.pi  # degrees version
        self.meshgrid['hv_meshd'] = self.meshgrid['hv_mesh'] * 180 / math.pi  # degrees version

    def _build_cylindrical(self, core_location_tuple=None):
        """
        Converts cartesian coordinate attributes into cylindrical attributes, and

        :param core_location_tuple:   tuple (X mm, Y mm) of actual core location on the meshgrid
        """

        # makes an initial guess at the core location, then builds the cylindrical coordinate system.
        if core_location_tuple is None:
            self._get_cylindrical_meshgrids(self._guess_core())

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
        for i in range(0, self.dims[-1]):  # cant figure out fully vectorized element wise subtraction
            self.dynamic_set['r'][:, :, i] = self.dynamic_set['R'][:, :, i] - self.mean_set['R']
            self.dynamic_set['t'][:, :, i] = self.dynamic_set['T'][:, :, i] - self.mean_set['T']

        self.mean_set['r'] = masked_rms(self.dynamic_set['r'], axis=2, mask=mpm)
        self.mean_set['t'] = masked_rms(self.dynamic_set['t'], axis=2, mask=mpm)

        # find dynamic reynolds stresses and turbulence
        for component in ['rr', 'tt', 'rt', 'rw', 'tw']:
            self.dynamic_set[component] = self.dynamic_set[component[0]] * self.dynamic_set[component[1]]
            self.mean_set[component] = masked_rms(self.dynamic_set[component], axis=2, mask=mpm)
        return

    def _get_circulation_strength(self, vtheta_max, core_radius):
        """ simply solves for the circulation strength of the vortex using Ash-Zardkahan model"""
        gamma = vtheta_max * 4 * math.pi * core_radius
        self.circulation_strength = gamma
        return gamma

    def get_spatial_derivatives_cylindrical(self):
        """
        Function obtains spatial derivatives in cylindrical coordinates. Using the chain
        rule to convert spatial derivatives from dfdx, dfdy into dfdr, and dfdt.
        :return:
        """

        # use the divide by zero function (dbz) to chain rule out the partial derivatives
        # shorten notation for r and t coordinates and convert to meters from mm
        r = self['r_mesh'] / 1000
        t = self['t_mesh']

        # chain rule the radial velocity derivatives
        drdx, drdy = get_spatial_derivative(self['R'], self['x_mesh'], self['y_mesh'])
        self.derivative_set['drdr'] = (dbz(drdx, np.cos(t)) +
                                       dbz(drdy, np.sin(t))) / 2
        self.derivative_set['drdt'] = (dbz(drdx, (1 / r) * np.cos(t)) +
                                       dbz(drdy, (-1 / r) * np.sin(t))) / 2

        # chain rule the tangential velocity derivatives
        dtdx, dtdy = get_spatial_derivative(self['T'], self['x_mesh'], self['y_mesh'])
        self.derivative_set['dtdr'] = (dbz(dtdx, np.cos(t)) +
                                       dbz(dtdy, np.sin(t))) / 2
        self.derivative_set['dtdt'] = (dbz(dtdx, (1 / r) * np.cos(t)) +
                                       dbz(dtdy, (-1 / r) * np.sin(t))) / 2

        # chain rule the axial velocity derivatives
        dwdx, dwdy = get_spatial_derivative(self['W'], self['x_mesh'], self['y_mesh'])
        self.derivative_set['dwdr'] = (dbz(dwdx, np.cos(t)) +
                                       dbz(dwdy, np.sin(t))) / 2
        self.derivative_set['dwdt'] = (dbz(dwdx, (1 / r) * np.cos(t)) +
                                       dbz(dtdy, (-1 / r) * np.sin(t))) / 2

        # assume indeterminite z partial derivatives are zero.
        self.derivative_set['drdz'] = self['W'] * 0
        self.derivative_set['dtdz'] = self['W'] * 0
        self.derivative_set['dwdz'] = self['W'] * 0

        # correction function!
        # there needs to be a correction function here!

        return self.derivative_set

    def get_x_axis_strain_rates(self):
        """
        This function obtains strain rate terms in CARTESIAN coordinates that approximate
        cylindrical coordinates on the x axis because here. r ~ x and y ~ theta.
        :return:
        """

        r = self['r_mesh'] / 1000
        x = self['x_mesh'] / 1000
        y = self['y_mesh'] / 1000
        t = self['t_mesh']

        # get the first term
        a = r * self['T']
        dadx, _ = get_spatial_derivative(a, x, y)
        term1, _ = get_spatial_derivative(dadx / r, x, y)
        self['sr_cx1'] = term1
        print 'term1', term1.max(), term1.min()

        # second term
        _, dtdy = get_spatial_derivative(self['T'], x, y * r)
        _, dtdy2 = get_spatial_derivative(dtdy, x, y * r)
        term2 = dtdy2 / (r * r)
        self['sr_cx2'] = term2

        # third term (unknown)
        term3 = term2 * 0
        self['sr_cx3'] = term3

        # fourth term
        _, drdy = get_spatial_derivative(self['R'], x, y)
        term4 = 2 * drdy / (r * r)
        self['sr_cx4'] = term4

        self.derivative_set['sr_cx_tot'] = term1 + term2 + term3 + term4
        return self.derivative_set

    def get_cylindrical_strain_rates(self):

        r = self['r_mesh'] / 1000
        x = self['x_mesh'] / 1000
        y = self['y_mesh'] / 1000
        t = self['t_mesh']

        # get the first term
        a = r * self['T']
        dadx, dady = get_spatial_derivative(a, x, y)
        dadr = (dbz(dadx, np.cos(t)) + dbz(dady, np.sin(t))) / 2
        dadr /= r
        dbdx, dbdy = get_spatial_derivative(dadr, x, y)
        term1 = (dbz(dbdx, np.cos(t)) + dbz(dbdy, np.sin(t))) / 2

        self['sr_cx1'] = term1
        print 'term1', term1.max(), term1.min()


        return self.derivative_set

    def get_pressure_relax_terms(self, pressure_relaxation=None):
        """
        Calculates the relatinshib between reynolds stress and turbulent viscosity as derived by the
        pressure relaxation equations.

        :param pressure_relaxation: the pressure relaxation coefficientin microseconds
        :return:
        """

        if pressure_relaxation is None:
            pressure_relaxation = 1

        # top = 1/r*2 * d/dr(r^2 * vtheta_vr reynolds stress)
        # bottom = second derivative of vtheta with respect to r +  d/dr(vtheta / r)
        # result = top / bottom

        r = self['r_mesh'] / 1000
        x = self['x_mesh'] / 1000
        y = self['y_mesh'] / 1000
        t = self['t_mesh']

        # this section of code with proper converted cylindrical derivatives has lots of opportunities
        # to divide by zero, so the results are ridiculously exploded in every case.
        def radius_chain_rule(quantity, x, y, r, t):
            """ takes radial derivatives by using coordinate transform chain rule """
            dqdx, dqdy = get_spatial_derivative(quantity, x, y)
            dqdr = (dbz(dqdx, np.cos(t)) + dbz(dqdy, np.sin(t))) / 2
            return dqdr

        top_dadr = radius_chain_rule(self['rt'], x, y, r, t)
        reynolds_term = top_dadr / (r ** 2)

        bot_dTdr = radius_chain_rule(self['T'], x, y, r, t)
        bot_dTdr2 = radius_chain_rule(bot_dTdr, x, y, r, t)
        bot_dTrdr = radius_chain_rule(self['T'] / r, x, y, r, t)
        vel_grad_term = bot_dTdr2 + bot_dTrdr

        drrdr = radius_chain_rule(self['rr'], x, y, r, t)
        noneq_press_term = self['T'] * (self['rr'] - self['tt'] + drrdr) / (r * r)
        self.equation_terms['turb_visc_reynolds'] = reynolds_term
        self.equation_terms['turb_visc_vel_grad'] = vel_grad_term
        self.equation_terms['turb_visc_ettap'] = (pressure_relaxation / 1.0e6) * noneq_press_term
        self.equation_terms['turb_visc_total'] = abs(dbz(noneq_press_term + reynolds_term, vel_grad_term))

        # ratio of the noneq turb visc and classical turb visc, no real relationship appears to be there
        self.equation_terms['turb_visc_ratio'] = (self.equation_terms['turb_visc_total'] /
                                                  self.equation_terms['turb_visc'])

        # now calculate the pressure gradient
        # this is closely related to the velocity gradient term
        mu = AIR_DYNAMIC_VISCOSITY
        dPdr = - (mu * vel_grad_term) / (pressure_relaxation * self['T'] / r)
        self.equation_terms['dPdr'] = dPdr

        '''
        # this section of code uses an approximate where x ~ r and  r * y ~ theta
        drsdr, _ = get_spatial_derivative(r * r * self['rt'], r, r * t)
        top = drsdr / (r ** 2)

        dTdr1, _ = get_spatial_derivative(self['T'], r, r * t)
        dTdr2, _ = get_spatial_derivative(dTdr1, r, r * t)
        dTrdr, _ = get_spatial_derivative(self['T'] / r, r, r * t)
        bot = dTdr2 + dTrdr

        self.equation_terms['turb_visc_reynolds'] = top
        self.equation_terms['turb_visc_vel_grad'] = bot
        self.equation_terms['turb_visc_ettap'] = dbz(top, bot)
        '''
        return


# ============== plotting functions=======================
    def get_smoothed_line(self, x, y, numpoints=None, M=None, std=None, order=None,
                          r_range=None, t_range=None, symmetric=None):
        """
        :param x:               either a component name string or a flattened array
        :param y:               either a component name string or a flattened array
        :param numpoints:       length of vectors, proportional to resolution
        :param M:               number of points in gaussian resample
        :param std:             standard deviation of points in gaussian resample
        :param order:           the number of times to repeat the smoothing process.
        """

        if numpoints is None:
            numpoints = 1e4
        if M is None:
            M = 20
        if std is None:
            std = 50
        if order is None:
            order = 100

        if isinstance(x, str):
            x = self._getitem_by_rt(x, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()
        if isinstance(y, str):
            y = self._getitem_by_rt(y, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()

        # smooth the curve
        x_avg, y_avg = smooth_filt(x, y, numpoints, M, std, convolve_mode='same', order=order)

        # chop off the end where strange things tend to happen.
        x_avg = x_avg[0:-int(len(x_avg) / 20)]
        y_avg = y_avg[0:-int(len(y_avg) / 20)]
        return x_avg, y_avg

    def _get_cbar_levels(self, component):
        comp = self._getitem_by_rt(component).astype('float')
        return np.nanpercentile(comp.filled(np.nan), range(1, 100))

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

    def _get_vrange(self, component, low_percentile=None, high_percentile=None,
                    r_range=None, t_range=None, symmetric=None):
        """
        Gets the percentile range for color bar scaling on a given component matrix. Used
        to ensure spurious high and low values do not over stretch the color ramps on plots.
        The more noise in the data the further form 0 and 100 respectively these percentiles
        must be.

        :param component:           component string to get range for, or a numpy array.
        :param low_percentile:      low percentile value marking coolest color
        :param high_percentile:     high percentile value marking warmest color
        :return:
        """

        if low_percentile is None:
            low_percentile = VRANGE_DEFAULT[0]
        if high_percentile is None:
            high_percentile = VRANGE_DEFAULT[1]

        if isinstance(component, str):
            comp = self._getitem_by_rt(component, r_range=r_range, t_range=t_range, symmetric=symmetric).astype('float')
        else:
            comp = component
        vmin = np.nanpercentile(comp.filled(np.nan), low_percentile)
        vmax = np.nanpercentile(comp.filled(np.nan), high_percentile)
        return vmin, vmax

    @staticmethod
    def _save_or_show(outpath=None):
        """
        :param outpath: output filepath to save figure, if left None, figure will be displayed to screen
        """
        if outpath is not None:
            plt.savefig(outpath, format='jpg', dpi=DEFAULT_DPI)
            print("saved figure to {0}".format(outpath))
            plt.close()
        else:
            plt.show()

    def _draw_core(self, fig, ax, normalized=False, color=None):
        """ draws a vortex core for reference on any contour plot """

        # default color of white
        if color is None:
            color = DEFAULT_CORE_MARKER_COLOR

        # plot the core location for reference
        if normalized:
            if self.core_location[0] is not None:
                ax.scatter(0, 0, marker='+', s=100, c=color)

            if self.core_radius is not None:
                circ = plt.Circle((0, 0), radius=1, edgecolor=color,
                                  linestyle=':', facecolor='none', label="Core Boundary")
        else:
            if self.core_location[0] is not None:
                ax.scatter(*self.core_location, marker='+', s=100, c=color)

            if self.core_radius is not None:
                circ = plt.Circle(self.core_location, radius=self.core_radius, edgecolor=color,
                                  linestyle=':', facecolor='none', label="Core Boundary")
        ax.add_patch(circ)

    def get_dvt_dr(self, t_range=None, r_range=None, symmetric=None, outpath=None):
        """
        Finds the derivative of Vtheta with respect to radius, a noteworthy quantity by sampling the
        space of T and r, creating a list ordered by r, applying a moving average, then finding the
        derivative.
        :return:
        """

        # pull in flattened arrays then sort them in ascending order
        r = self._getitem_by_rt('r_mesh', r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()
        t = self._getitem_by_rt('T', r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()

        # take moving averages for smoothing before taking derivative
        ravg, tavg = smooth_filt(r, t, 1e4, 100, 50, order=50)

        # chop of back 5% to remove the smoothing distortion
        ravg=ravg[0:-int(len(ravg) / 20)]
        tavg=tavg[0:-int(len(tavg) / 20)]

        rdr, dtavgdr = get_spatial_derivative(tavg, ravg)
        dtavgdr *= 1000  # convert from mm to m

        # find the maximum point along the moving average fit for vortex characterization
        t_max = tavg.max()
        r_t_max = ravg[np.unravel_index(np.ma.argmax(tavg), tavg.shape)]

        # set up plots with sizing and labels
        fig, ax1 = plt.subplots(figsize=(8, 4), dpi=DEFAULT_DPI, facecolor='w')
        smooth_label = shorthand_to_tex('T')
        deriv_label = "$\\frac{{\\partial \\overline{t}}}{{\partial R}}$"

        # plot the first axis in terms of tangential velocity
        ax1.scatter(r / r_t_max, t, marker='.', color='black', s=0.5)
        smooth = ax1.plot(ravg / r_t_max, tavg, linestyle='-', color='blue', label=smooth_label)
        ax1.set_ylabel(shorthand_to_tex('T'))
        ax1.set_xlabel("Core Radii")
        ax1.set_ylim(0)

        # on a second axis plot the derivative curve in green
        ax2 = ax1.twinx()
        deriv = ax2.plot(rdr / r_t_max, dtavgdr, linestyle='--', color='green', label=deriv_label)
        ax2.set_ylabel(deriv_label)
        plt.grid()
        plt.legend()

        plt.xlim(0, 4)
        plt.tight_layout()

        if outpath is not False:
            self._save_or_show(outpath)
        else:
            plt.close()

        results_dict = {"r_t_max": r_t_max,
                        "t_max": t_max,
                        "r": r,
                        "t": t,
                        "ravg": ravg,
                        "tavg": tavg}
        return results_dict

    def comparison_plot(self, pressure_relaxation=None, t_range=None, r_range=None,
                        symmetric=None, outpath=None, kinematic_viscosity=None, time=None):
        """
        This function is very similar to get_dvt_dr, and is just another step in our evolving understanding
        of what we actually want to show here. It plots the average experimental profile against several
        theoretical vortex profiles including Rankin, Lamb-Oseen, and Ash.

        :param pressure_relaxation:     the pressure relaxation coefficient in microseconds
        :param t_range:                 theta range to fit scatter data
        :param r_range:                 r range to fit scatter data
        :param symmetric:               use symmetry about each axis
        :param outpath:                 filepath to save the plot as an image
        :param kinematic_viscosity:     the kinematic viscosity in units of m^2 / s (read from config file)
        :param time:                    the age of the vortex in seconds, used for Lamb-Oseen decay simulation
        :return:
        """

        if kinematic_viscosity is None:
            kinematic_viscosity = AIR_KINEMATIC_VISCOSITY

        # pull in flattened arrays then sort them in ascending order (convert to meters)
        r_scatter = self._getitem_by_rt('r_mesh', r_range=r_range,
                                        t_range=t_range, symmetric=symmetric).flatten() / 1000
        t_scatter = self._getitem_by_rt('T', r_range=r_range,
                                        t_range=t_range, symmetric=symmetric).flatten()

        # take moving averages for smoothing before taking derivative
        r_array, t_exp = smooth_filt(r_scatter, t_scatter, 1e4, 100, 50, order=50)

        # chop of back 5% to remove the smoothing distortion
        r_array = r_array[0:-int(len(r_array) / 20)]
        t_exp = t_exp[0:-int(len(t_exp) / 20)]

        # find the maximum point along the moving average fit for vortex characterization
        vtheta_max = t_exp.max()
        core_radius = r_array[np.unravel_index(np.ma.argmax(t_exp), t_exp.shape)]

        # get the circulation strength
        circulation_strength = self._get_circulation_strength(vtheta_max, core_radius)

        # rankine vortex
        t_rankine = RankineVortex(core_radius, circulation_strength).get_vtheta(r_array)

        # lamb oseen vortex based on core_radius
        lamboseen = LambOseenVortex(circulation_strength, kinematic_viscosity)
        t_lamboseen = lamboseen.get_vtheta(r_array, vtheta_max=vtheta_max, core_radius=core_radius)

        # ash vortex based on vtheta max
        t_ash = AshVortex(core_radius, vtheta_max=vtheta_max).get_vtheta(r_array)

        # burnham and hallock vortex
        #t_burnhamhallock = BurnhamHallockVortex(core_radius, circulation_strength).get_vtheta(r_array)

        # set up plots with sizing and labels
        fig, ax1 = plt.subplots(figsize=(8, 4), dpi=DEFAULT_DPI, facecolor='w')

        r_plot = r_array / core_radius
        r_scat = r_scatter / core_radius
        plt.scatter(r_scat, t_scatter, marker='.', color='lightgray', s=0.3, label="Experimental Data")
        plt.plot(r_plot, t_exp, linewidth=2, linestyle='-', color='black', label="Experimental Fit")
        plt.plot(r_plot, t_rankine, linewidth=2, linestyle='--', color='firebrick', label="Rankine: $r_{core}$")
        plt.plot(r_plot, t_lamboseen, linewidth=2, linestyle=':', color='navy', label="Lamb-Oseen: $V_{\\theta, max}$")
        plt.plot(r_plot, t_ash, linewidth=2, linestyle='-.', color='purple', label="Ash et al.: $V_{\\theta, max}$")

        if pressure_relaxation is not None:     # ash vortex based on pressure relaxation
            ettap  = pressure_relaxation / 1e6
            t_ash2 = AshVortex(core_radius, circulation_strength, kinematic_viscosity, ettap).get_vtheta(r_array)
            plt.plot(r_plot, t_ash2, linewidth=2, linestyle='-.', color='darkgreen', label="Ash: $\etta_p$")

        plt.legend()
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)
        plt.xlabel("$R/r_{core}$")
        plt.ylabel(shorthand_to_tex('T'))
        plt.title("Theoretical azimuthal velocity profile fits", fontsize=DEFAULT_TITLE_SIZE)

        if outpath is not False:
            self._save_or_show(outpath)
        else:
            plt.close()
        return

    def dynamic_plot(self, component_y, title=None, y_label=None, cmap=None,
                     r_range=None, t_range=None, symmetric=None, tight=False,
                     figsize=None, outpath=None):
        """
        Creates multi plot showing time history of dynamic data over all runs, plus spectral content.

        :param component_y:     the component under study
        :param title:           title for the plot
        :param y_label:         plots y axis label
        :param cmap:            colormap to use for scatter plot
        :param r_range:         the range of radius to sample for dynamic averaging
        :param t_range:         the tangential range to sample for dynamic averaging
        :param symmetric:       set "True" if radial and tangential components should be symetric about xy
        :param tight:           set to True for a tightly ploted area
        :param figsize:         figure size in inches at 200 dpi
        :param outpath:         output path to save the plot as a eps.
        :return:
        """

        r_range = self._rrange_parser(r_range)

        if component_y not in DYNAMIC_INCLUDES:
            raise Exception("That component is not saved in the dynamic dataset! check DYNAMIC_INCLUDES")

        print("generating dynamic plot...")
        if title is None:
            title = "{0} vs time".format(shorthand_to_tex(component_y))
        if y_label is None:
            y_label = shorthand_to_tex(component_y)
        if figsize is None:
            figsize = (10, 5)
        if cmap is None:
            cmap = SCATTER_DEFAULT_CMAP
        x_label = "Time (seconds)"

        subset_kwargs = {'r_range': r_range, 't_range': t_range, 'symmetric': symmetric}
        y_sets, t_set = self._get_dynamic_subsets(component_y, **subset_kwargs)

        # now make the figure
        fig = plt.figure(figsize=figsize, dpi=DEFAULT_DPI, facecolor='w')

        # first plot
        gs = plt.GridSpec(100, 100, bottom=0.15, left=0.02, right=0.98)
        ax1 = fig.add_subplot(gs[:, 5:24])

        xplot_mesh = (self['x_mesh'] - self.core_location[0]) / self.core_radius
        yplot_mesh = (self['y_mesh'] - self.core_location[1]) / self.core_radius
        areaplot = ax1.contourf(xplot_mesh, yplot_mesh,
                                self._getitem_by_rt(component_y, **subset_kwargs),
                                CONTOUR_DEFAULT_LEVELS,
                                cmap=CONTOUR_DEFAULT_CMAP)

        # create the colorbar with just a few axis ticks. nticks should just be a kwarg for colorbar, idk why it isn't.
        ticks = list(self._get_vrange(component_y, r_range=r_range, t_range=t_range, symmetric=symmetric,
                                      low_percentile=0, high_percentile=100))
        cb = plt.colorbar(areaplot, orientation="horizontal", ticks=ticks, pad=0.2)

        plt.gca().set_aspect('equal', adjustable='box')         # set equal aspect ratio, x/y
        self._draw_core(None, ax1, normalized=True, color='k')  # add the circular core marker.

        plt.title("Sample Area", fontsize=DEFAULT_TITLE_SIZE - 2)
        plt.xlabel("$X/r_{core}$")
        plt.ylabel("$Y/r_{core}$")
        plt.xlim(-r_range[1] / self.core_radius, r_range[1] / self.core_radius)
        plt.ylim(-r_range[1] / self.core_radius, r_range[1] / self.core_radius)
        plt.locator_params(nbins=3)     # controls number of x and y axis ticks to avoid crowding

        # second dynamic plot with percentile bounds as well.
        ax2 = fig.add_subplot(gs[:, 30:77])
        ax2.plot(t_set, y_sets['mean'], 'k-', label='Mean')
        ax2.plot(t_set, y_sets['p05'], 'k:', label='90% of Values')
        ax2.plot(t_set, y_sets['p95'], 'k:')
        plt.title(title, fontsize=DEFAULT_TITLE_SIZE - 2)
        plt.grid()
        plt.legend(loc=1)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # third plot psd
        ax3 = fig.add_subplot(gs[:, 84:100])
        ax3.psd(y_sets['mean'], NFFT=64, Fs=SAMPLING_RATE)
        plt.title("log PSD", fontsize=DEFAULT_TITLE_SIZE - 2)

        self._save_or_show(outpath)
        return {"t_set": t_set, "y_sets": y_sets}

    def scatter_plot(self, component_x, component_y, component_c=None,
                     title=None, x_label=None, y_label=None, c_label=None, cmap=None,
                     x_range=None, y_range=None, r_range=None, t_range=None, symmetric=None,
                     tight=True, figsize=None, outpath=None, log_y=None, show_grid=False):
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
        :param outpath:         filepath to save the figure
        :param log_y:           set True to make the y axis logarithmic
        :param show_grid:       set True to show dotted grid on the plot
        :param smooth:          set True to smooth the data, may also math smooth_kwargs
        :param smooth_kwargs:   keyword arguments to pass to the smoothing function, if any. see get_average_profile
        :return:
        """

        # gather argument values which are left as defaults.
        if title is None:
            title = "{0} vs {1}".format(shorthand_to_tex(component_y), shorthand_to_tex(component_x))
        if x_label is None:
            x_label = shorthand_to_tex(component_x)
        if y_label is None:
            y_label = shorthand_to_tex(component_y)
            if log_y:
                y_label = "Log of {0}".format(y_label)
        if c_label is None and component_c is not None:
            c_label = shorthand_to_tex(component_c)
        if figsize is None:
            figsize = (8, 4)
        if cmap is None:
            cmap = SCATTER_DEFAULT_CMAP

        # get the scatter datasets.
        x = self._getitem_by_rt(component_x, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()
        y = self._getitem_by_rt(component_y, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()

        # convert the data to logarithmic
        if log_y is None:
            log_y = False
        elif log_y is True:
            y = abs(y)

        # if the x axis is the core radius, normalize by the radius of the core.
        if component_x == "r_mesh":
            x /= self.core_radius
            x_label = "$R/r_{core}$"

        fig = plt.figure(figsize=figsize, dpi=DEFAULT_DPI, facecolor='w', edgecolor='k')
        if component_c is not None:
            c = self._getitem_by_rt(component_c, r_range=r_range, t_range=t_range, symmetric=symmetric).flatten()
            vmin, vmax = self._get_vrange(component_c, r_range=r_range, t_range=t_range)
            plt.scatter(x, y, marker=SCATTER_DEFAULT_MARKER, c=c, cmap=cmap, vmax=vmax, vmin=vmin)
            cb = plt.colorbar(orientation='vertical')
            cb.set_label(c_label)
        else:
            plt.scatter(x, y, marker=SCATTER_DEFAULT_MARKER, color=SCATTER_DEFAULT_COLOR)
            c = None

        # if smoothing is desired, get smoothed versions of the datasets
        # not added

        # apply manual specifications of x and y range, otherwise guess.
        if y_range is None:
            vmin, vmax = self._get_vrange(y, 0, 100)
            plt.ylim(vmin - 0.1, vmax * 1.1)
        else:
            plt.ylim(y_range[0], y_range[1])
        if x_range is not None:
            plt.xlim(x_range[0], x_range[1])
        if tight:
            plt.tight_layout()
        if show_grid:
            plt.grid()
        if log_y:
            plt.yscale('log')

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title, fontsize=DEFAULT_TITLE_SIZE, y=1.04)

        self._save_or_show(outpath)

        # return the plotted data
        return {"x": x, "y": y, "c": c}

    def scatter_plot_qual(self, component_x, component_y):
        """
        Prints quick simple scatter plot of component_x vs component_y with the points colored
        according to the number of samples making up data from that point. Useful for evaluating
        trends and differentiating between real trends and potentially spurious features.

        params are exactly as scatter_plot()
        """
        return self.scatter_plot(component_x, component_y, 'num', c_label="Number of Samples")

    def quiver_plot(self, title=None, outpath=None):
        """
        creates a quiver plot of the vector field
        """

        if title is None:
            title = "Quiver plot"

        plt.figure()
        plt.quiver(self['x_mesh'], self['y_mesh'],
                   self['U'], self['V'],
                   color=SCATTER_DEFAULT_COLOR,
                   scale=400,
                   width=0.001,
                   headwidth=2,
                   headlength=1,
                   minshaft=1)

        plt.title(title, fontsize=DEFAULT_TITLE_SIZE)
        plt.tight_layout()
        xlims, ylims = self._get_plot_lims(50, 50)
        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.xlabel("$X$ (mm)")
        plt.ylabel("$Y$ (mm)")
        self._save_or_show(outpath)

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
        plt.title(title, fontsize=DEFAULT_TITLE_SIZE)

        # plot the core location for reference
        if self.core_location[0] is not None:
            ax.scatter(*self.core_location, marker='+', s=200, c='black')

        plt.tight_layout()
        xlims, ylims = self._get_plot_lims()
        plt.xlabel("$X/r_{core}$")
        plt.ylabel("$Y/r_{core}$")
        plt.xlim(xlims)
        plt.ylim(ylims)
        plt.axis('equal')

        self._save_or_show(outpath)
        return

    def contour_plot(self, component, t_range=None, r_range=None, symmetric=False,
                     title=None, outpath=None, log_colorbar=False, diverging=False):
        """
        creates a contour plot of input component

        :param component:   component to plot, any member of the meshgrid or the mean_set
        :param title:       custom title to place on the contour plot
        :param outpath:     an output path to save a eps file of this plot
        :return:
        """

        # sets the color ramp to diverge
        if diverging:
            cmap = CONTOUR_DIVERGE_CMAP
        else:   # standard color ramp
            cmap = CONTOUR_DEFAULT_CMAP

        if isinstance(component, np.ma.MaskedArray):
            data = component
        else:
            data = self._getitem_by_rt(component, r_range=r_range, t_range=t_range, symmetric=symmetric)

        if title is None:
            title = shorthand_to_tex(component)

        xplot_mesh = (self['x_mesh'] - self.core_location[0]) / self.core_radius
        yplot_mesh = (self['y_mesh'] - self.core_location[1]) / self.core_radius
        vmin, vmax = self._get_vrange(component, r_range=r_range, t_range=t_range, symmetric=symmetric)
        fig, ax = plt.subplots(figsize=(8, 7), dpi=DEFAULT_DPI)
        if log_colorbar:
            cf = plt.contourf(xplot_mesh, yplot_mesh, data, CONTOUR_DEFAULT_LEVELS,
                              norm=LogNorm(), cmap=cmap, vmin=vmin, vmax=vmax)
        else:
            cf = plt.contourf(xplot_mesh, yplot_mesh, data, CONTOUR_DEFAULT_LEVELS,
                              cmap=cmap, vmin=vmin, vmax=vmax)

        cf.set_clim(vmin=vmin, vmax=vmax)
        plt.colorbar(cf)
        plt.title(title, fontsize=DEFAULT_TITLE_SIZE, y=1.04)
        plt.xlabel("$X/r_{core}$")
        plt.ylabel("$Y/r_{core}$")
        plt.xlim(-5, 5)
        plt.ylim(-5, 5)

        self._draw_core(fig, ax, normalized=True)
        plt.tight_layout()
        plt.gca().set_aspect('equal', adjustable='box')

        self._save_or_show(outpath)
        return

    def hist_plot(self, component, bins=100):
        """ plots quick and dirty histogram of a component """
        c = self._getitem_by_rt(component).flatten()
        fig = plt.figure()
        plt.hist(c[~c.mask], bins=bins, color="grey")
        plt.show()
        return




if __name__ == "__main__":

    DEFAULT_DPI = 120
    exp_num = 55
    force = False
    if not os.path.exists("temp{0}.pkl".format(exp_num)) or force:
        directory = os.path.join(DATA_FULL_DIR, str(exp_num))
        paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".v3d")]
        mvf = AxialVortex("temp{0}".format(exp_num), v3d_paths=paths, min_points=20)
        mvf.find_core()
        mvf.get_cart_turbulent_viscosity()
        mvf.get_pressure_relax_terms()
        mvf.to_pickle("temp{0}.pkl".format(exp_num), include_dynamic=False)
    else:
        mvf = AxialVortex().from_pickle("temp{0}.pkl".format(exp_num))
        mvf.get_pressure_relax_terms()




    #mvf.scatter_plot('r_mesh', 'turb_visc') #, log_y=True)

    kwargs = {"t_range": (10, 80),
              "r_range": (0, '3r'),
              "symmetric": True,
              "show_grid": True,
              #"y_range": (1e1, 1e9),
              }

    from matplotlib import cm

    #mvf.contour_plot('dPdr', r_range=('1r','5r'), t_range=(20,70), symmetric=True, cmap=cm.PRGn)
    #mvf.contour_plot('dPdr', r_range=('0.3r','3r'), t_range=(20,70), symmetric=True, cmap=cm.PRGn)
    mvf.scatter_plot('r_mesh', 'dPdr', 'ctke', r_range=('1r','5r'), t_range=(20,70), symmetric=True)
    #mvf.scatter_plot('r_mesh', 'turb_visc_reynolds', log_y=True, **kwargs)
    #mvf.scatter_plot('r_mesh', 'turb_visc_vel_grad', log_y=True, **kwargs)
    #mvf.scatter_plot('r_mesh', 'turb_visc_ettap', log_y=True, **kwargs)
    #mvf.scatter_plot('r_mesh', 'turb_visc_total', **kwargs)
    #mvf.scatter_plot('r_mesh', 'turb_visc_ratio')


