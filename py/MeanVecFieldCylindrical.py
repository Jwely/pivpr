__author__ = 'Jwely'

from MeanVecField import MeanVecField


class MeanVecFieldCylindrical(MeanVecField):
    def __init__(self, v3d_paths, name_tag):
        # invoke the parent class init
        MeanVecField.__init__(self, v3d_paths, name_tag)

        # add cylindrical specific attributes
        self.core_location = (None, None)       # position of core
        self.core_index = (None, None)          # fractional index position of core

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

        # all the same data in the matrices above, but flattened into a 1d array and in terms of "distance to center"
        self.vel_flat = {'U': None,  # x direction mean velocity
                         'V': None,  # y direction mean velocity
                         'W': None,  # z direction mean velocity
                         'R': None,  # mean radial velocity around vortex core
                         'T': None,  # mean tangential velocity around vortex core

                         'u': None,  # U fluctuation
                         'v': None,  # V fluctuation
                         'w': None,  # W fluctuation
                         'r': None,  # fluctuation in R
                         't': None,  # fluctuation in T

                         'uu': None,  # turbulent energy in u
                         'vv': None,  # turbulent energy in v
                         'ww': None,  # turbulent energy in w
                         'cte': None,  # total cartesian turbulent energies
                         'rr': None,  # turbulent energy in radial
                         'tt': None,  # turbulent energy in tangential

                         'uv': None,  # reynolds stress in u/v
                         'uw': None,  # reynolds stress in u/w
                         'vw': None,  # reynolds stress in v/w
                         'crs': None,  # total cartesian reynolds stresses
                         'rt': None,  # reynolds stress in r/t
                         'rw': None,  # reynolds stress in r/w
                         'tw': None,  # reynolds stress in t/w
                         'num': None}  # number of good data points making up stats for other values


if __name__ == "__main__":
    paths = [r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01000.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01001.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01002.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01003.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01004.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01005.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01006.v3d",
             r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01007.v3d"]

    mvf = MeanVecFieldCylindrical(paths, "test")
    mvf.show('num')