__author__ = 'Jwely'

""" little module for defining global variables """

from matplotlib import cm
from matplotlib import rcParams

# this prevents plt.tight_layout() from crowding axis labels off the edges of the plot.
rcParams.update({'figure.autolayout': True})
print("configured tight layout!")

# global variables
CONTOUR_DEFAULT_LEVELS = 32         # number of discreet colors to use on contour plots
CONTOUR_DEFAULT_CMAP = cm.Greys     # default color ramp to use on contour plots
CONTOUR_DEFAULT_RRANGE = (0, 100)   # default radius range to subset contour plots by

SCATTER_DEFAULT_CMAP = cm.jet       # default colormap to use on scatter plots with third variable
SCATTER_DEFAULT_RLIM = 100          # default plot



