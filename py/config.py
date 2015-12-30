__author__ = 'Jwely'

""" little module for defining global variables """

from matplotlib import cm
from matplotlib import rcParams

# this prevents plt.tight_layout() from crowding axis labels off the edges of the plot.
rcParams.update({'figure.autolayout': True})
print("configured tight layout!")

# global variables
CONTOUR_DEFAULT_LEVELS = 16
CONTOUR_DEFAULT_CMAP = cm.Greys
CONTOUR_DEFAULT_RRANGE = (0, 100)
SCATTER_DEFAULT_CMAP = cm.jet


