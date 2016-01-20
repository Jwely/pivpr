__author__ = 'Jwely'

from os.path import join, dirname, abspath
from matplotlib import cm
from matplotlib import rcParams
config_root = dirname(abspath(__file__))


# this prevents plt.tight_layout() from crowding axis labels off the edges of the plot.
rcParams.update({'figure.autolayout': True})


# resource filepaths and directories
EXPERIMENT_TABLE_PATH = abspath(join(config_root, "piv/dat/experiment_table.csv"))
PIV_PICKLE_DIR = abspath(join(config_root, "piv/pickles"))
DATA_FULL_DIR = abspath(join(config_root, "../data_full"))
CALIBRATION_DIR = abspath(join(config_root, "uncertainty/cal_data"))
SYNTHESIZED_PIV_DIR = abspath(join(config_root, "uncertainty/artificial_images"))

TEX_FIGURE_DIR = abspath(join(config_root, "../texdocs/figs"))
TEX_TABLE_DIR = abspath(join(config_root, "../texdocs/tables"))
TEX_MAIN_PATH = abspath(join(config_root, "../texdocs/main.tex"))


# when pickling, save only these dynamic components (saves disk space), discard the others
DYNAMIC_INCLUDES = ['ctke', 'r', 't', 'w', 'rt', 'rw', 'tw']


# statistics variables
DEFAULT_MIN_POINTS = 20             # minimum number of points required to consider a data point as good


# global plotting variables
CONTOUR_DEFAULT_LEVELS = 32         # number of discreet colors to use on contour plots
CONTOUR_DEFAULT_CMAP = cm.Greys     # default color ramp to use on contour plots (greyscale friendly)
CONTOUR_DEFAULT_RRANGE = (0, 100)   # default radius range to subset contour plots by
SCATTER_DEFAULT_CMAP = cm.Greys     # default colormap to use on scatter plots with third variable
SCATTER_DEFAULT_RLIM = 100          # default plot
VRANGE_DEFAULT = (1, 99.9)          # default percentile values for defining color ramp boundaries

# attributes of the dataset
SAMPLING_RATE = 1                   # the sampling rate in Hz. (1/time between vector fields)



