__author__ = "Jwely"

import math
import numpy


class LambOseenVortex:
    def __init__(self, circulation_strength, viscosity):
        """
        Creates a Rankine vortex with points across
        :param circulation_strength:    The circulation strength of the vortex
        :param viscosity:               The kinematic viscosity in units of (m^2/s)
                                        is about 0.154 for air at room temperature
        :return:
        """

        self.circulation_strength = circulation_strength
        self.viscosity = viscosity
        pass


    def get_vtheta(self, r_array, core_radius=None, time=None, verbose=True):
        """
        returns array of vtheta values for input r_array and other parameters. uses
        core_radius input if available, otherwise a time input is required.

        :param r_array:     array of radial points for which to compute vtheta
        :param core_radius: the radius of the core in meters!
        :param time:        the age of the vortex in seconds, or downstream distance / free stream
        :param verbose:     set False to suppress print statements.
        :return:
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

        if verbose:
            print("Lamb-Oseen vortex core radius at time {0}s is {1}m".format(time, core_radius))

        return vtheta