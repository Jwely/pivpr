__author__ = 'Jwely'

from py.constructors import build_axial_vortex


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
    av.get_contour_plot('T')


if __name__ == "__main__":
    test_data()
