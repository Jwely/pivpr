__author__ = 'Jwely'

""" little module for defining global variables """

from matplotlib import cm
from matplotlib import rcParams

# this prevents plt.tight_layout() from crowding axis labels off the edges of the plot.
rcParams.update({'figure.autolayout': True})
print("configured tight layout!")

# when pickling, save only these dynamic components (saves disk space), discard the others
DYNAMIC_INCLUDES = ['ctke', 'r', 't', 'w', 'rt', 'rw', 'tw']

# statistics variables
DEFAULT_MIN_POINTS = 20             # minimum number of points required to consider a data point as good

# global plotting variables
CONTOUR_DEFAULT_LEVELS = 32         # number of discreet colors to use on contour plots
CONTOUR_DEFAULT_CMAP = cm.jet       # default color ramp to use on contour plots (greyscale friendly)
CONTOUR_DEFAULT_RRANGE = (0, 100)   # default radius range to subset contour plots by
SCATTER_DEFAULT_CMAP = cm.jet       # default colormap to use on scatter plots with third variable
SCATTER_DEFAULT_RLIM = 100          # default plot
VRANGE_DEFAULT = (1, 99.9)          # default percentile values for defining color ramp boundaries

# attributes of the dataset
SAMPLING_RATE = 1                   # the sampling rate in Hz. (1/time between vector fields)



