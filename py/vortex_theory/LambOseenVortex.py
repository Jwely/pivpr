__author__ = "Jwely"

import math
import numpy


class LambOseenVortex:
    def __init__(self, circulation_strength, viscosity):
        """
        Creates a Lamb Oseen vortex
        :param circulation_strength:    The circulation strength of the vortex
        :param viscosity:               The kinematic viscosity in units of (m^2/s)
                                        is about 0.154 for air at room temperature
        :return:
        """

        self.circulation_strength = circulation_strength
        self.viscosity = viscosity
        pass

    def _get_vtheta_by_vthetamax(self, r_array, vtheta_max, core_radius):
        """
        Alternative method to obtain lamb-oseen profile by Vtheta_max
        :param vtheta_max:  the maximum azimuthal velocity
        :param core_radius: the radius of the core in meters!
        """

        a = 1.25643
        rmax = a ** 0.5 * core_radius

        vtheta = vtheta_max * (1 + 1 / (2 * a)) * (rmax / r_array) * (1 - numpy.exp(-a * r_array ** 2 / rmax ** 2))
        return vtheta


    def _get_vtheta_standard(self, r_array, core_radius=None, time=None):
        """
        returns array of vtheta values for input r_array and other parameters. uses
        core_radius input if available, otherwise a time input is required.

        :param r_array:     array of radial points for which to compute vtheta
        :param core_radius: the radius of the core in meters!
        :param time:        the age of the vortex in seconds, or downstream distance / free stream
        """

        gamma = self.circulation_strength

        if core_radius is not None:
            time = core_radius / (4 * self.viscosity)       # compute time
        elif time is not None:
            core_radius = (4 * self.viscosity * time)       # compute core radius
        else:
            raise Exception("must use core_radius or time parameter")

        # find vtheta
        exponent = (-1.0 * (r_array / core_radius) ** 2)
        vtheta = gamma / (2 * math.pi * r_array) * (1 - numpy.exp(exponent))
        return vtheta


    def get_vtheta(self, r_array, vtheta_max=None, core_radius=None, time=None, verbose=True):
        """
        Get the vtheta profile with one of a few input combinations
        1) with vtheta_max and core_radius
        2) with core radius only
        3) with the age of the vortex

        the first option relies upon nothing except the inputs, the second and third options depend
        upon the time variant nature of the vortex and the viscosity and circulation strength attributes
        of the LambOseenVortex instance.

        :param r_array:     array of radial points for which to compute vtheta
        :param vtheta_max:  the maximum azimuthal velocity
        :param core_radius: the radius of the core in meters!
        :param time:        the age of the vortex in seconds, or downstream distance / free stream
        :param verbose:     set to false to suppress information outputs.
        :return:
        """
        if vtheta_max is None:
            vtheta = self._get_vtheta_standard(r_array, core_radius, time)
        elif vtheta_max is not None and core_radius is not None:
            vtheta = self._get_vtheta_by_vthetamax(r_array, vtheta_max, core_radius)
        else:
            raise Exception("Invalid input combination")

        if verbose:
            print("Lamb-Oseen vortex core radius at time {0}s is {1}m".format(time, core_radius))

        return vtheta