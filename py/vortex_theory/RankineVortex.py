__author__ = "Jwely"

import math


class RankineVortex:

    def __init__(self, core_radius, circulation_strength):
        """
        Creates a Rankine vortex
        :param core_radius: the core radius in (meters)
        :param circulation_strength: the circulation strength of the vortex
        :return:
        """

        self.circulation_strength = circulation_strength
        self.core_radius = core_radius
        pass


    def get_vtheta(self, r_array):
        """
        Returns an array of vtheta values for input r_array
        :param r_array: array of radius points at which to sample vtheta
        :return:
        """

        gamma = self.circulation_strength
        rcore = self.core_radius

        # find vtheta on inside and outside of core boundary
        vtheta_in = gamma / (2 * math.pi * rcore) * (r_array / rcore)
        vtheta_out = gamma / (2 * math.pi * r_array)

        # assemble piecewise output vtheta
        vtheta = vtheta_out
        vtheta[r_array < rcore] = vtheta_in[r_array < rcore]

        return vtheta
