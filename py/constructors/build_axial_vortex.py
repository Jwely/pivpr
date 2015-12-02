__author__ = 'Jwely'

import os
from py.managers import AxialVortex


def build_axial_vortex(v3d_dir, pkl_dir, name_tag, include_dynamic=False,
                       velocity_fs=None, force_recalc=False):
    """
    Returns an AxialVortex instance. Manages pickling of classes to allow one-time-computation
    as often as possible, unless the user overrides it with force_recalc=True.

    :param v3d_dir:         directory containing v3d files containing data on axial vortex
    :param pkl_dir:         directory containing binary pickle files
    :param name_tag:        a unique name identifier for the AxialVortex instance (used in pkl filename)
    :param include_dynamic: set to True to keep matrix data from every individual run instead of taking averages
                            and standard deviations then throwing out the individual datasets. Substantially
                            increases system requirements to keep True.
    :param velocity_fs:     The value of the free stream velocity in meters/second
    :param force_recalc:    set True to force re-computation of all parameters from raw datasets.
    :return:
    """

    # generate pkl_path based on the pkl_dir and name_tag inputs.
    if include_dynamic:
        pkl_path = os.path.join(os.path.abspath(pkl_dir), "{0}_dyn.pkl".format(name_tag))
    else:
        pkl_path = os.path.join(os.path.abspath(pkl_dir), "{0}.pkl".format(name_tag))


    if os.path.exists(pkl_path) and not force_recalc:
        av_instance = AxialVortex().from_pickle(pkl_path)
        return av_instance

    else:
        v3d_paths = [os.path.join(v3d_dir, fname) for fname in os.listdir(v3d_dir) if fname.endswith(".v3d")]
        av_instance = AxialVortex(name_tag=name_tag, v3d_paths=v3d_paths, velocity_fs=velocity_fs)
        av_instance.find_core()                           # find core location (not yet proven 100% robust)
        av_instance.build_cylindrical()                   # convert cartesian data to cylindrical about core
        av_instance.to_pickle(pkl_path, include_dynamic)  # pickle it so we can just load it again later
        return av_instance


if __name__ == "__main__":

    kwargs = {"v3d_dir": r"E:\Data2\Ely_May28th\Vector\1",
              "pkl_dir": "../pickles",
              "name_tag": "01",
              "include_dynamic": True,
              "velocity_fs": 15.22,
              "force_recalc": False}

    av = build_axial_vortex(**kwargs)
    av.show_contour('P')