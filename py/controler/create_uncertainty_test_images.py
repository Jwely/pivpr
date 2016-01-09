__author__ = 'Jwely'

import os
from py.uncertainty import ArtificialPIV


def create_uncertainty_test_images():
    """
    Builds a bunch of test images for monte carlo uncertainty analysis.


    NOTE: This takes a while. The process of computing light intensities from particles
    and adding them to the final output image is probably very inefficient, but I don't know
    how to make it better and can't be bothered to figure it out. I'd rather just hit run and let
    it go overnight.

    :return:
    """

    uncertainty_image_dir = "../uncertainty/artificial_images"

    velocity_sets = {'high_z': {'dt': 25, 'x_vel': 0.00, 'y_vel': 0.00, 'z_vel': 40.0},  # high z
                     'low_z': {'dt': 40, 'x_vel': 0.00, 'y_vel': 0.00, 'z_vel': 10.0},  # low z
                     'high_x': {'dt': 25, 'x_vel': 10.0, 'y_vel': 0.00, 'z_vel': 0.00},  # high x
                     'mid_x': {'dt': 25, 'x_vel': 1.00, 'y_vel': 0.00, 'z_vel': 0.00},  # mid x
                     'low_x': {'dt': 40, 'x_vel': 0.10, 'y_vel': 0.00, 'z_vel': 0.00},  # low x
                     'high_y': {'dt': 25, 'x_vel': 0.00, 'y_vel': 10.0, 'z_vel': 0.00},  # high y
                     'mid_y': {'dt': 25, 'x_vel': 0.00, 'y_vel': 1.00, 'z_vel': 0.00},  # mid y
                     'low_y': {'dt': 40, 'x_vel': 0.00, 'y_vel': 0.10, 'z_vel': 0.00}}  # low y

    calibrations = [1, 7]

    dimensions = (1024, 1280)

    for vel_key in velocity_sets.keys():
        for cal in calibrations:
            name = "loc-{0}-{1}".format(cal, vel_key)

            outdir = os.path.join(uncertainty_image_dir, name)
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            cal_path = "../uncertainty/cal_data/station_{num}/ely_may28th.cal".format(num=cal)
            apiv = ArtificialPIV(dimensions, name)
            apiv.load_calibration_file(cal_path)
            apiv.make_image_pairs(n_particles=dimensions[0] * dimensions[1] / 40,
                                  output_dir=outdir,
                                  particle_size=0.2,
                                  particle_scatter=100,
                                  light_sheet_thickness=3.0,
                                  **velocity_sets[vel_key])


if __name__ == "__main__":
    create_uncertainty_test_images()
