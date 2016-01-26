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


    def get_vtheta(self, r_array, time, verbose=True):
        """
        Returns an array of vtheta values for input r_array and time
        :param r_array:     Array of radius points at which to sample vtheta
        :param time:        Age of the vortex in seconds. For axial vortices, this is the
                            downstream distance divided by the free stream velocity
        :return:
        """

        gamma = self.circulation_strength
        rcore = (4 * self.viscosity * time)     # compute core radius as function of time and viscosity

        # find vtheta
        exponent = (- (r_array ** 2) / rcore ** time)
        vtheta = gamma / (2 * math.pi * r_array) * (1 - numpy.exp(exponent))

        if verbose:
            print("Lamb-Oseen vortex core radius at time {0}s is {1}m".format(time, rcore))

        return vtheta