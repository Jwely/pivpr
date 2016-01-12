__author__ = 'Jwely'

import os
from py.uncertainty import ArtificialPIV
from py.config import *
import copy


def synthesize_piv_uncertainty_images(uncertainty_image_dir):
    """
    Builds a bunch of test images for monte carlo uncertainty analysis.


    NOTE: This takes a while. The process of computing light intensities from particles
    and adding them to the final output image is probably very inefficient, but I don't know
    how to make it better and can't be bothered to figure it out. I'd rather just hit run and let
    it go over night (scratch that ... several days).
    """
    '''
    velocity_sets = {'high_z': {'dt': 25, 'u': 0.00, 'v': 0.00, 'w': 40.0},  # high z
                     'low_z': {'dt': 40, 'u': 0.00, 'v': 0.00, 'w': 10.0},  # low z
                     'high_x': {'dt': 25, 'u': 10.0, 'v': 0.00, 'w': 0.00},  # high x
                     'mid_x': {'dt': 25, 'u': 1.00, 'v': 0.00, 'w': 0.00},  # mid x
                     'low_x': {'dt': 40, 'u': 0.10, 'v': 0.00, 'w': 0.00},  # low x
                     'high_y': {'dt': 25, 'u': 0.00, 'v': 10.0, 'w': 0.00},  # high y
                     'mid_y': {'dt': 25, 'u': 0.00, 'v': 1.00, 'w': 0.00},  # mid y
                     'low_y': {'dt': 40, 'u': 0.00, 'v': 0.10, 'w': 0.00}}  # low y
    '''
    velocity_sets = {'high25': {'dt': 25, 'u': 10.0, 'v': 10.0, 'w': 40.0},
                     'low25': {'dt': 25, 'u': 0.05, 'v': 0.05, 'w': 15.0},
                     'high40': {'dt': 40, 'u': 6.0, 'v': 6.0, 'w': 28.0},
                     'low40': {'dt': 40, 'u': 0.05, 'v': 0.05, 'w': 5.00}}

    calibrations = [1, 2, 3, 4, 5, 6, 7]

    dimensions = (1024/2, 1280/2)

    for vel_key in velocity_sets.keys():
        for cal in calibrations:
            name = "loc-{0}-{1}".format(cal, vel_key)

            outdir = os.path.join(uncertainty_image_dir, name)
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            cal_path = os.path.join(CALIBRATION_DIR, "station_{num}/ely_may28th.cal".format(num=cal))
            apiv = ArtificialPIV(dimensions, name)
            apiv.load_calibration_file(cal_path)

            velocity_sets[vel_key].update({"n_particles": dimensions[0] * dimensions[1] / 40,
                                           "output_dir": outdir,
                                           "particle_size": 0.2,
                                           "particle_scatter": 100,
                                           "light_sheet_thickness": 3.0})

            apiv.make_image_pairs(**velocity_sets[vel_key])




if __name__ == "__main__":
    synthesize_piv_uncertainty_images("../uncertainty/artificial_images")
