__author__ = 'Jwely'


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


class VecField3d:


    def __init__(self, filepath):
        """
        All meaningful attributes of a VecField3d object are constructed upon __init__.
        The "output" attributes is the velocity_matrix
        :param filepath:        local filepath to a .v3d file
        """

        self.filepath = filepath                # local filepath to .v3d file
        self.headers = None                     # list of column headers
        self.dataframe = None                   # pandas dataframe of csv like data.

        # set up empty coordinate dictionary
        self.dims = (None, None)                # x, y dimensions of all matrix data

        # set up empty coordinate value dictionary, (x and y are two-dimensionalized 1d vectors)
        self.x_set = None
        self.y_set = None
        self.meshgrid = {"x": None,
                         "y": None}

        self.velocity_matrix = {'U': None,      # x direction velocity
                                'V': None,      # y direction velocity
                                'W': None}      # z direction velocity

        # Build up the attributes
        self._read_v3d()                        # parse the .v3d file, populate self.headers, self.dataframe
        self._get_meshgrid()                    # populate self.dims, self.meshgrid
        self._table_to_matrix()                 # populate self.velocity_matrix
        print("loaded {0}".format(filepath))


    def _read_v3d(self):
        """ parses the .v3d file into a pandas dataframe """

        # get the header row
        with open(self.filepath, "r") as f:
            header_raw = f.next()

        # parse the header_raw string to extract list of header strings
        self.headers = header_raw.split("VARIABLES=")[1].replace('"', "").replace("\n", "").split(", ")

        # now load the rest of the file with a pandas dataframe
        self.dataframe = pd.read_csv(self.filepath, skiprows=1, names=self.headers)


    def _get_meshgrid(self):
        """
        Fills the following attributes of the object:
            self.meshgrid
            self.dims
        """

        self.x_set = sorted(set(self.dataframe['X mm']))
        self.y_set = sorted(set(self.dataframe['Y mm']))
        self.dims = (len(self.y_set), len(self.x_set))

        y_mesh, x_mesh = np.meshgrid(self.y_set, self.x_set)
        self.meshgrid['x'] = x_mesh
        self.meshgrid['y'] = y_mesh

        for component in ['U', 'V', 'W']:
            self.velocity_matrix[component] = np.zeros(self.dims)


    def _table_to_matrix(self):
        """
        Fills matrix values for cartesian mesh. These values tend to be organized
        """

        for i, row in self.dataframe.iterrows():
            x_index = self.x_set.index(row['X mm'])
            y_index = self.y_set.index(row['Y mm'])

            self.velocity_matrix['U'][y_index, x_index] = row['U m/s']
            self.velocity_matrix['V'][y_index, x_index] = row['V m/s']
            self.velocity_matrix['W'][y_index, x_index] = row['W m/s']

        # turn the arrays into masked arrays to hide nodata
        for component in ['U', 'V', 'W']:
            masked = np.ma.masked_array(self.velocity_matrix[component],
                                        mask=self.velocity_matrix[component] > 100)
            self.velocity_matrix[component] = masked


    def show(self):

        for component in ['U', 'V', 'W']:
            plot_data = self.velocity_matrix[component]
            fig, ax = plt.subplots()
            heatmap = ax.pcolor(plot_data)
            plt.show()



if __name__ == "__main__":

    fpath = r"E:\Data2\Ely_May28th\Vector\1\Ely_May28th01000.v3d"
    v = VecField3d(fpath)
    v.show()

