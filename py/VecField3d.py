__author__ = 'Jwely'


import pandas as pd
import numpy as np



class VecField3d:

    def __init__(self, filepath):
        """
        All meaningful attributes of a VecField3d object are constructed upon __init__.
        The "output" attributes is the velocity_matrix
        :param filepath:    local filepath to a .v3d file
        """
        self.filepath = filepath                # local filepath to .v3d file
        self.headers = None                     # list of column headers
        self.dataframe = None                   # pandas dataframe of csv like data.

        # set up empty coordinate dictionary
        self.dims = (None, None)                # x, y dimensions of all matrix data

    # each of the attributes below are actually dictionaries, where the key is the dimension and the
    # value is a matrix with dimensions equal to self.dims

        # set up empty coordinate value dictionary, (x and y are two-dimensionalized 1d vectors)
        self.cords_matrix = {"x": None,
                             "y": None}

        self.velocity_matrix = {"u": None,      # x direction
                                "v": None,      # y direction
                                "w": None,      # z (stream-wise) direction
                                "r": None,      # radial about vortex core
                                "t": None}      # tangential about vortex core

        self.reynolds_matrix = {'uu': None,     # turbulent energy in u
                                'vv': None,     # turbulent energy in v
                                'ww': None,     # turbulent energy in w
                                'uv': None,     # reynolds stress in u/v
                                'uw': None,     # reynolds stress in u/w
                                'vw': None,     # reynolds stress in v/w

                                'rr': None,     # turbulent energy in radial
                                'tt': None,     # turbulent energy in tangential
                                'rt': None,     # reynolds stress in r/t
                                'rw': None,     # reynolds stress in r/w
                                'tw': None,     # reynolds stress in t/w
                                'tot': None}    # total turbulent energies

        # all the same data in the matrices above, but flattened into a 1d array and in terms of "distance to core"
        self.velocity_flat = {"u": None,        # x direction
                              "v": None,        # y direction
                              "w": None,        # z (streamwise) direction
                              "r": None,        # radial about vortex core
                              "t": None}        # tangential about vortex core

        self.reynolds_flat = {'uu': None,       # turbulent energy in u
                              'vv': None,       # turbulent energy in v
                              'ww': None,       # turbulent energy in w
                              'uv': None,       # reynolds stress in u/v
                              'uw': None,       # reynolds stress in u/w
                              'vw': None,       # reynolds stress in v/w

                              'rr': None,       # turbulent energy in radial
                              'tt': None,       # turbulent energy in tangential
                              'rt': None,       # reynolds stress in r/t
                              'rw': None,       # reynolds stress in r/w
                              'tw': None,       # reynolds stress in t/w
                              'tot': None}      # total turbulent energies

    # Build up the attributes


        self._read_v3d()                        # parse the .v3d file
        self._get_matrix_dims()                 # populate self.dims, coord_matrix





    def _read_v3d(self):
        """ parses the .v3d file into a pandas dataframe """

        # get the header row
        with open(self.filepath, "r") as f:
            header_raw = f.next()

        # parse the header_raw string to extract list of header strings
        self.headers = header_raw.split("VARIABLES=")[1].replace('"', "").replace("\n", "").split(", ")

        # now load the rest of the file with a pandas dataframe
        self.dataframe = pd.read_csv(self.filepath, skiprows=1, names=self.headers)


    def _get_matrix_dims(self):

        self.x_set = sorted(set(self.dataframe['X mm']))
        self.y_set = sorted(set(self.dataframe['Y mm']))

        print(len(self.x_set))
        print(len(self.y_set))




if __name__ == "__main__":

    fpath = r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01000.v3d"
    v = VecField3d(fpath)


