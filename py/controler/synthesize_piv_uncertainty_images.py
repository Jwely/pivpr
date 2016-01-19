__author__ = 'Jwely'

import os
from py.uncertainty import ArtificialPIV
from py.config import *
import copy


def synthesize_piv_uncertainty_images(uncertainty_image_dir):
    """
    Builds a bunch of test images for monte carlo uncertainty analysis. This function
    creates outputs specifically to work with Insight PIV software by TSI.

    NOTE: This takes a while. The process of computing light intensities from particles
    and adding them to the final output image is probably very inefficient, but I don't know
    how to make it better and can't be bothered to figure it out. I'd rather just hit run and let
    it go over night (scratch that ... several days).
    """

    velocity_order = ['min', 'low', 'high', 'max']
    velocity_sets = {'min': {'dt': 40, 'u': 0.0, 'v': 0.0, 'w': 5.0},
                     'low': {'dt': 40, 'u': 4.0, 'v': 4.0, 'w': 26.0},
                     'high': {'dt': 25, 'u': 0.0, 'v': 0.0, 'w': 9.0},
                     'max': {'dt': 25, 'u': 10.0, 'v': 10.0, 'w': 38.0}
                     }

    calibrations = [7]

    dimensions = (1024, 1280)

    for i, vel_key in enumerate(velocity_order):
        for cal in calibrations:
            name = "Ely_May28th{0}{1}".format(str(cal).zfill(2), str(i + 1).zfill(3))

            cal_path = os.path.join(CALIBRATION_DIR, "station_{num}/ely_may28th.cal".format(num=cal))
            apiv = ArtificialPIV(dimensions, name)
            apiv.load_calibration_file(cal_path)

            velocity_sets[vel_key].update({"n_particles": dimensions[0] * dimensions[1] / 15,
                                           "particle_size": 0.2,
                                           "particle_scatter": 100,
                                           "light_sheet_thickness": 3.0,
                                           "dtype": "uint16",
                                           "name": name,
                                           "output_dir": uncertainty_image_dir})

            apiv.make_image_pairs(**velocity_sets[vel_key])



if __name__ == "__main__":
    synthesize_piv_uncertainty_images("../uncertainty/artificial_images")
    import datetime
    print datetime.datetime.now()
