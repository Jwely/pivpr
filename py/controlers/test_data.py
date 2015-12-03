__author__ = 'Jwely'

from py.constructors import build_axial_vortex

# some of my favorite colormaps for quick reference
from matplotlib import cm
DIV_CMAP = cm.seismic           # diverging colormap
SEQ_CMAP = cm.bone_r            # sequential colormap
CYC_CMAP = cm.hsv               # cyclical colormap
DEF_CMAP = cm.jet               # default colormap
QUA_CPAP = cm.nipy_spectral     # qualitative spectrum colormap


def test_data():

    kwargs = {"v3d_dir": "../../test_data",
              "pkl_dir": "../pickles",
              "name_tag": "test_data",
              "include_dynamic": False,
              "velocity_fs": 15.22,
              "force_recalc": False}

    av = build_axial_vortex(**kwargs)
    #av.get_stream_plot()
    #av.get_contour_plot('P')
    av.get_scatter_plot2('r_mesh', 'W')
    av.get_scatter_plot('r_mesh', 'T', 'cte', cmap=cm.jet,
                        title="Tangential Velocity and Turbulent Kinetic Energy",
                        x_label="Distance to vortex core (mm)",
                        y_label="Tangential velocity (m/s)",
                        c_label="Turbulent kinetic energy (2TKE)")
    #av.get_scatter_plot2('r_mesh', 'cte')
    #av.get_contour_plot('T')


if __name__ == "__main__":
    test_data()
