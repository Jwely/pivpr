__author__ = 'Jwely'

from os.path import join, dirname, abspath
from matplotlib import cm
from matplotlib import rcParams
config_root = dirname(abspath(__file__))


# this prevents plt.tight_layout() from crowding axis labels off the edges of the plot.
rcParams.update({'figure.autolayout': True})

# aerodynamic params to assume
AIR_DENSITY = 1.225                     # kg/m^3
AIR_KINEMATIC_VISCOSITY = 15.68e-6      # m^2 / s
AIR_DYNAMIC_VISCOSITY = 18.46e-6        # kg / m*s

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
DYNAMIC_INCLUDES = ['U', 'V', 'ctke', 'r', 't', 'w', 'rt', 'rw', 'tw']


# statistics variables
DEFAULT_MIN_POINTS = 15             # minimum number of points required to consider a data point as good


# global plotting variables
CONTOUR_DEFAULT_LEVELS = 256        # number of discreet colors to use on contour plots
CONTOUR_DEFAULT_CMAP = cm.jet       # default color ramp to use on contour plots
CONTOUR_DIVERGE_CMAP = cm.PRGn      # default color ramp for diverging contour plots (things centered about zero)
CONTOUR_DEFAULT_RRANGE = (0, 50)    # default radius range to subset contour plots by
SCATTER_DEFAULT_CMAP = cm.jet       # default colormap to use on scatter plots with third variable
SCATTER_DEFAULT_COLOR = 'navy'      # default color used on scatter plots with just one variable
DEFAULT_CORE_MARKER_COLOR = 'k'     # default color of lines and crosshair used to mark the core boundary
SCATTER_DEFAULT_MARKER = 'x'        # marker to use on scatter plots
SCATTER_DEFAULT_RLIM = 100          # default plot
VRANGE_DEFAULT = (1, 99)            # default percentile values for defining color ramp boundaries
SCATTER_VRANGE_DFAULT = (5, 95)      # default percentile value for defining scatter plot y axis
DEFAULT_DPI = 300                   # default dots per inch
DEFAULT_TITLE_SIZE = 18             # font size for titles

# attributes of the dataset
SAMPLING_RATE = 1                   # the sampling rate in Hz. (1/time between vector fields)


