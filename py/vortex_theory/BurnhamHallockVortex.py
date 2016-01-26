__author__ = "Jwely"

import math


class BurnhamHallockVortex:
    def __init__(self, core_radius, circulation_strength):
        self.core_radius = core_radius
        self.circulation_strength = circulation_strength


    def get_vtheta(self, r_array):
        """
        Returns an array of vtheta values for input r_array. the method for this
        varies depending upon how the class was instantiated
        :param r_array:  array of radius values for which vtheta values are wanted.
        """

        gamma = self.circulation_strength
        rcore = self.core_radius
        vtheta = gamma / (2 * math.pi * r_array) * (r_array ** 2 / (r_array ** 2 + rcore ** 2))
        return vtheta





