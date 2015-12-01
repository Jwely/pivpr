__author__ = 'Jwely'

from MeanVecFieldCartesian import MeanVecFieldCartesian
import numpy as np


class MeanVecFieldCylindrical(MeanVecFieldCartesian):

    def __init__(self, v3d_paths, name_tag):
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

        :param v3d_paths:       list of filepaths to v3d files
        :param name_tag:        unique string name tag for this data set
        :return:
        """

        # invoke the parent class init
        MeanVecFieldCartesian.__init__(self, v3d_paths, name_tag)

        # add cylindrical specific attributes
        self.core_location = (None, None)       # position of core
        self.core_index = (None, None)          # fractional index position of core

        self.meshgrid.update({"r_mesh": None,   # radial meshgrid
                              "t_mesh": None})  # tangential meshgrid

        # add to the velocity matrix and flattened version
        self.vel_matrix.update({'R': None,  # mean radial velocity around vortex core
                                'T': None,  # mean tangential velocity around vortex core

                                'r': None,  # fluctuation in R
                                't': None,  # fluctuation in T

                                'rr': None,  # turbulent energy in radial
                                'tt': None,  # turbulent energy in tangential

                                'rt': None,  # reynolds stress in r/t
                                'rw': None,  # reynolds stress in r/w
                                'tw': None})  # reynolds stress in t/w


    def define_core(self, index_location=None, real_location=None):
        """
        Accepts a manual definition of x,y coordinates of the core location.

        :param index_location:  tuple of (x,y) index locations (may be floats for fractional index positioning)
        :param real_location:   tuple of (X mm, Y mm) real coordinate locations
        :return:
        """

        pass


    def _cartesian_to_cylindrical(self, real_location):
        """
        Converts cartesian coordinate attributes into cylindrical attributes, and

        :param real_location:   tuple (X mm, Y mm) of actual core location on the meshgrid
        :return:
        """

        xc, yc = real_location

        # r = ((x - x_core)**2 + (y - y_core)**2)**0.5
        # t = numpy.math.atan2((x - x_core), (y - y_core))

        # build the cylindrical meshgrid
        self['r_mesh'] = ((self['x_mesh'] - xc) ** 2 + (self['y_mesh'] - yc) ** 2) ** 0.5
        self['t_mesh'] = np.math.atan2((self['y_mesh'] - yc), (self['x_mesh'] - xc))


        # need to fill these vel_matrix attributes with this function
        blah = {'R': None,  # mean radial velocity around vortex core
                'T': None,  # mean tangential velocity around vortex core

                'r': None,  # fluctuation in R
                't': None,  # fluctuation in T

                'rr': None,  # turbulent energy in radial
                'tt': None,  # turbulent energy in tangential

                'rt': None,  # reynolds stress in r/t
                'rw': None,  # reynolds stress in r/w
                'tw': None}  # reynolds stress in t/w


    def show_scatter(self, component_y, component_x):
        """
        prints quick simple scatter plot of component_x vs component_y. Useful for viewing data
        as a function of distance to vortex core (R) or angle around the core (T) """

        # will need to 1dimensionalize the vel_matrix data for plotting purposes.
        pass


if __name__ == "__main__":
    paths = [r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01000.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01001.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01002.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01003.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01004.v3d",
             ]

    mvf = MeanVecFieldCylindrical(paths, "test")
    mvf.show_heatmap('num')
    mvf.define_core(index_location=(74, 61))