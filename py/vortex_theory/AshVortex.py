__author__ = "Jwely"

import math


class AshVortex:
    def __init__(self, core_radius, circulation_strength=None, viscosity=None,
                 pressure_relaxation=None, vtheta_max=None):
        """
        Based on work by Robert Ash, Irfan Zardadkhan, Allan Zuckerwar in
        "The influence of pressure relaxation on the structure of an axial vortex"

        inputs must either be:
            core_radius
            circulation_strength
            viscosity
            pressure_relaxation
        OR
            core_radius
            vtheta_max

        :param core_radius:
        :param circulation_strength:
        :param viscosity:
        :param pressure_relaxation:
        :param vtheta_max:
        :return:
        """

        self.circulation_strength = circulation_strength
        self.pressure_relaxation = pressure_relaxation
        self.viscosity = viscosity
        self.vtheta_max = vtheta_max
        self.core_radius = core_radius

        if vtheta_max is not None:
            self.method = "vtheta_max"
        elif circulation_strength is not None and viscosity is not None and pressure_relaxation is not None:
            self.method = "circulation_strength"
        else:
            raise Exception("incomplete input combination, read the docstring!")


    def _get_vtheta_other(self, r_array):
        """ handles unsimplified form """

        gamma = self.circulation_strength
        rcore = self.core_radius
        ettap = self.pressure_relaxation
        nu = self.viscosity

        Rgamma = gamma / (2 * math.pi * nu)
        vtheta = (gamma / math.pi) * (2.0 ** 1.5) / (Rgamma * (nu * ettap) ** 0.5) * \
                 (r_array / rcore) / ((r_array / rcore) ** 2 + 1)
        return vtheta


    def _get_vtheta_from_max(self, r_array):
        """ handles simplified form in which vtheta max is already known """

        vtmax = self.vtheta_max
        rcore = self.core_radius

        vtheta = 2 * vtmax * (r_array / rcore) / ((r_array / rcore) ** 2 + 1)
        return vtheta


    def get_vtheta(self, r_array, verbose=True):
        """
        Returns an array of vtheta values for input r_array. the method for this
        varies depending upon how the class was instantiated
        """

        if verbose:
            if self.method == "vtheta_max":
                print("Solving with user supplied vtheta_max")
            elif self.method == "circulation_strength":
                print("Solving with circulation strength, viscosity, and pressure relaxation inputs")

        if self.method == "vtheta_max":
            return self._get_vtheta_from_max(r_array)
        elif self.method == "circulation_strength":
            return self._get_vtheta_other(r_array)




