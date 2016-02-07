__author__ = 'Jwely'

import os
from py.piv import AxialVortex


def construct_axial_vortex(v3d_dir, pkl_dir, name_tag, include_dynamic=False,
                           velocity_fs=None, z_location=None, min_points=20, force_recalc=False):
    """
    Returns an AxialVortex instance. Manages pickling of classes to allow one-time-computation
    as often as possible, unless the user overrides it with force_recalc=True.

    :param v3d_dir:         directory containing v3d files containing data on axial vortex
    :param pkl_dir:         directory containing binary pickle files
    :param name_tag:        a unique name identifier for the AxialVortex instance (used in pkl filename)
    :param include_dynamic: set to True to keep matrix data from every individual run instead of taking
                            averages and standard deviations then throwing out the individual datasets.
                            Substantially increases system requirements to keep True.
    :param velocity_fs:     The value of the free stream velocity in meters/second
    :param z_location:      the location downstream from the vortex generator in mm
    :param force_recalc:    set True to force re-computation of all parameters from raw datasets.

    :return AxialVortex:
    """

    # generate pkl_path based on the pkl_dir and name_tag inputs.
    pkl_path_dyn = os.path.join(os.path.abspath(pkl_dir), "{0}_dyn.pkl".format(name_tag))
    pkl_path = os.path.join(os.path.abspath(pkl_dir), "{0}.pkl".format(name_tag))

    if include_dynamic:
        if os.path.exists(pkl_path_dyn) and not force_recalc:
            av_instance = AxialVortex().from_pickle(pkl_path_dyn)
            return av_instance
    else:
        if os.path.exists(pkl_path) and not force_recalc:
            av_instance = AxialVortex().from_pickle(pkl_path)
            return av_instance

    # build list of all v3d files in directory from which to create the AxialVortex
    v3d_paths = [os.path.join(v3d_dir, fname) for fname in os.listdir(v3d_dir) if fname.endswith(".v3d")]
    av_instance = AxialVortex(name_tag=name_tag, v3d_paths=v3d_paths,
                              velocity_fs=velocity_fs, z_location=z_location, min_points=min_points)

    # find the core and build cylindrical coordinate data around it
    av_instance.find_core()
    av_instance.calculate_turbulent_viscosity()

    # pickle both the full dynamic set and the reduced set
    av_instance.to_pickle(pkl_path_dyn, include_dynamic=True)
    av_instance.to_pickle(pkl_path, include_dynamic=False)
    return av_instance