# coding=utf-8
__author__ = 'Jwely'

import numpy as np
import os
import json
from datetime import datetime
from py.utils.tiff_tools import save_array_as_dtype


class Particle:
    def __init__(self, x_mm, y_mm, z_mm):
        self.x_mm = x_mm
        self.y_mm = y_mm
        self.z_mm = z_mm


class ArtificialPIV:
    """
    This class is used to create simulated PIV images with very precisely calculated
    particle displacements. Particle displacements are propagated forward onto the
    camera image planes with intensity functions.
    """

    def __init__(self, dims, name=None):
        """
        :param dims:    a touple of the dimensions, ex: (1280, 1024)
        :param name:    a name for the class instance, used for naming output files
        """

        if name is None:
            self.name = "unnamed"
        else:
            self.name = name  # a name for the instance, used in saving output files

        self.dims = dims  # dimension tuple of the instance
        self.particles_a = []  # a list that will be filled with particle instances
        self.particles_b = []  # a list that will be filled with particle instances

        x_mesh, y_mesh = np.meshgrid(range(dims[1]), range(dims[0]))

        self.mesh = {'x_px': x_mesh,            # integer meshgrid of x pixel indices
                     'y_px': y_mesh,            # integer meshgrid of y pixel indices
                     'z_mm': np.zeros(dims)}    # integer meshgrid of z indices (all zero)

        self.images = {'La': np.zeros(dims),  # left image at time = 0
                       'Lb': np.zeros(dims),  # left image at time = dt
                       'Ra': np.zeros(dims),  # right image at time = 0
                       'Rb': np.zeros(dims)}  # right image at time = dt

        # dictionary of actual calibration information
        self.cal_info = {'l_x_mm': [],  # from left pixels to interrogation plane mm
                         'l_y_mm': [],  # from left pixels to interrogation plane mm
                         'l_x_px': [],  # real displacement to left pixel displacement
                         'l_y_px': [],  # real displacement to left pixel displacement
                         'r_x_mm': [],  # from right pixels to interrogation plane mm
                         'r_y_mm': [],  # from right pixels to interrogation plane mm
                         'r_x_px': [],  # real displacement to right pixel displacement
                         'r_y_px': []}  # real displacement to right pixel displacement

        # simple dictionary to translate text in the uncertainty file into tags used for parsing
        self.tag_dict = {'Left Camera x mm': 'l_x_mm',
                         'Left Camera y mm': 'l_y_mm',
                         'Left Camera X pix': 'l_x_px',
                         'Left Camera Y pix': 'l_y_px',
                         'Right Camera x mm': 'r_x_mm',
                         'Right Camera y mm': 'r_y_mm',
                         'Right Camera X pix': 'r_x_px',
                         'Right Camera Y pix': 'r_y_px'}


    def load_calibration_file(self, filepath):
        """
        Parses a .uncertainty file of format 'TSI_CAL_VERSION 2.100000' and extracts the equation
        information. This is expressed as a list of 9 coefficients and the corresponding orders of the
        x,y,z displacement terms to be multiplied. So a row of [123, 0, 1, 1] would translate to '123 * y * z'

        :param filepath:
        :return:
        """

        def get_active(aline):
            # identify the active section for further parsing
            for key in self.tag_dict.keys():
                if key in aline:
                    return self.tag_dict[key]

            if aline.isalpha():
                return None
            return None

        active = None

        with open(filepath, 'r') as calfile:
            for line in calfile:

                if line[0] in ['L', 'R', 'C']:
                    active = get_active(line)
                    continue

                if active is not None:
                    coeff, xord, yord, zord = line.replace('\n', '').split(', ')
                    self.cal_info[active].append([float(coeff), int(xord), int(yord), int(zord)])
        return


    @staticmethod
    def _get_overlap_fov(mesh_dict):
        """ Finds the bounding rectangle of real coordinate space that is visible to both cameras """

        result = {'x_mm_min': max([np.min(mesh_dict['l_x_mm']), np.min(mesh_dict['r_x_mm'])]),
                  'x_mm_max': min([np.max(mesh_dict['l_x_mm']), np.max(mesh_dict['r_x_mm'])]),
                  'y_mm_min': max([np.min(mesh_dict['l_y_mm']), np.min(mesh_dict['r_y_mm'])]),
                  'y_mm_max': min([np.max(mesh_dict['l_y_mm']), np.max(mesh_dict['r_y_mm'])])}

        result.update({'x_mm_cen': np.mean([result['x_mm_min'], result['x_mm_max']]),
                       'y_mm_cen': np.mean([result['y_mm_min'], result['y_mm_max']])})
        return result


    def _eval_cal_equation(self, cal_param, x, y, z):
        """
        Simply evaluates a displacement value from the input x,y,z coordinates.
        Note that depending on the cal_param, x, y, z may have units of (mm) or (pix).

        This is not computationally efficient, but it's effective and concise.
        """

        if cal_param not in self.cal_info.keys():
            raise Exception('Invalid cal_param {0}'.format(cal_param))

        ev = 0.00
        for line in self.cal_info[cal_param]:
            coeff, xord, yord, zord = line
            ev += coeff * (x ** xord) * (y ** yord) * (z ** zord)
        return ev


    def get_pixel_coords(self, x_mm=None, y_mm=None, z_mm=None, particle=None):
        """
        gets pixel coordinates on the left and right cameras. Accepts either a set of
        x y and z coordinates or a Particle instance which has those attributes.
        :param x_mm:    x position in mm
        :param y_mm:    y position in mm
        :param z_mm:    z position in mm
        :param particle:    a Particle instance which has the other three params as attributes
        """

        if particle is not None:
            x_mm = particle.x_mm
            y_mm = particle.y_mm
            z_mm = particle.z_mm
            return {'l_x_px': self._eval_cal_equation('l_x_px', x_mm, y_mm, z_mm),
                    'l_y_px': self._eval_cal_equation('l_y_px', x_mm, y_mm, z_mm),
                    'r_x_px': self._eval_cal_equation('r_x_px', x_mm, y_mm, z_mm),
                    'r_y_px': self._eval_cal_equation('r_y_px', x_mm, y_mm, z_mm)}

        elif x_mm is not None and y_mm is not None and z_mm is not None:
            return {'l_x_px': self._eval_cal_equation('l_x_px', x_mm, y_mm, z_mm),
                    'l_y_px': self._eval_cal_equation('l_y_px', x_mm, y_mm, z_mm),
                    'r_x_px': self._eval_cal_equation('r_x_px', x_mm, y_mm, z_mm),
                    'r_y_px': self._eval_cal_equation('r_y_px', x_mm, y_mm, z_mm)}


    def get_mm_coords(self, x_px, y_px, z_mm):
        """
        gets real mm coordinates from pixel coordinates.
        :param x_px:    x position in pixels
        :param y_px:    y position in pixels
        :param z_mm:    z position in millimeters
        """
        result = {'l_x_mm': self._eval_cal_equation('l_x_mm', x_px, y_px, z_mm),
                  'l_y_mm': self._eval_cal_equation('l_y_mm', x_px, y_px, z_mm),
                  'r_x_mm': self._eval_cal_equation('r_x_mm', x_px, y_px, z_mm),
                  'r_y_mm': self._eval_cal_equation('r_y_mm', x_px, y_px, z_mm)}
        return result


    def _get_intensities(self, particles, mesh_dict, particle_size=0.2, particle_scatter=100,
                         light_sheet_thickness=3.0, subset_radius=2.0):
        """
        Calculate the intensity of light surrounding a particle. Based on work from Raffel et al.

        Raffel, M., Willert, C., and Kompenhans, J., 'Particle Image Velocimetry: A Practical
        Guide', Springer, New York, 1998.

        :param particles:               a list of particle instances with x, y, z mm positions
        :param mesh_dict:               a dictionary of meshgrids as output from self.get_mm_coords. This is very
                                        important to use, because without creating a new meshgrid for every
                                        z_mm particle displacement, variation in the z direction can be lost!
        :param particle_size:           the size of the typical particle (in pixels!)
        :param particle_scatter:        scattering efficiency of particle, between 0 and 100
        :param light_sheet_thickness:   thickness in mm of the laser light sheet. (between +3 / -3 sigma intensity)
        :param subset_radius:           the furthest radius (in mm) to consider for particle intensities.
        :return:                        left and right intensity matrices, to be added to the master.

        NOTE:
        About 70% of the processing time is eaten up by this function. Efficiency gains here are worth a lot.
        Some significant improvements have already been made. The conditionals used to subset the matrix into
        smaller pieces are writen for speed.
        """

        if not isinstance(particles, list):
            particles = [particles]

        i_l = np.zeros(mesh_dict['l_x_mm'].shape)
        i_r = np.zeros(mesh_dict['l_x_mm'].shape)

        for i, p in enumerate(particles):
            # intensity after particle scattering.
            i_o = particle_scatter * np.exp(-(p.z_mm ** 2) / (light_sheet_thickness ** 2) / 8)

            # performance enhancement. Only executes the gaussian intensity function within 'index_radius'
            # grid points. This is because the intensity trails of drastically with significant distance, and
            # calculating the intensity for every pixel in every image for every particle is many million
            # more operations than is reasonably needed. This moves the scale from O(n^2) towards O(n) if n is the
            # size of one side on a square image
            l_sub = abs(mesh_dict['l_x_mm'] - p.x_mm) < subset_radius
            l_sub *= abs(mesh_dict['l_y_mm'] - p.y_mm) < subset_radius

            r_sub = abs(mesh_dict['r_x_mm'] - p.x_mm) < subset_radius
            r_sub *= abs(mesh_dict['r_y_mm'] - p.y_mm) < subset_radius

            # intensities at each camera
            i_l[l_sub] += i_o * np.exp((-(mesh_dict['l_x_mm'][l_sub] - p.x_mm) ** 2 -
                                        (mesh_dict['l_y_mm'][l_sub] - p.y_mm) ** 2) /
                                       ((particle_size ** 2) / 8))
            i_r[r_sub] += i_o * np.exp((-(mesh_dict['r_x_mm'][r_sub] - p.x_mm) ** 2 -
                                        (mesh_dict['r_y_mm'][r_sub] - p.y_mm) ** 2) /
                                       ((particle_size ** 2) / 8))

            if i % (len(particles) / 20) == 0:
                print("{0:2.0f}% complete".format(float(i) / len(particles) * 100))

        return i_l, i_r


    def make_image_pairs(self, n_particles, dt, u, v, w,
                         particle_size=0.2, particle_scatter=100, light_sheet_thickness=3.00,
                         dtype="uint16", name=None, output_dir=None):
        """
        This function creates a set of image pairs at times 0 and dt. It does this by creating
        particles within view of both cameras, then propagating them forward to the image planes
        of each one. It then moves the same particles according to input velocity and time parameters
        and repeats the process. This generates a complete PIV data set with precise known velocity
        vectors which can be sent through PIV processing for uncertainty analysis.

        :param n_particles:             number of particles to spawn in the interrogation plane
        :param dt:                      time step between frames (microseconds)
        :param u:                       uniform x velocity to displace particles
        :param v:                       uniform y velocity to displace particles
        :param w:                       uniform z velocity to displace particles
        :param particle_size:           average size of particles in (mm)
        :param particle_scatter:        scattering efficiency of the particle
        :param light_sheet_thickness:   thickness of light sheet (+/-3 sigma) in (mm)
        :param dtype:                   datatype of output images, defaults to uint16
        :param name:                    this name will be used to name output files, followed by La,Lb,Ra,Rb
        :param output_dir:              directory to save output images
        :return:
        """
        start = datetime.now()

        if name is None:
            name = self.name

        # translate velocities into displacements
        x_disp = u * dt * 1e-6 * 1e3  # convert to mm from m/s and microseconds
        y_disp = v * dt * 1e-6 * 1e3  # convert to mm from m/s and microseconds
        z_disp = w * dt * 1e-6 * 1e3  # convert to mm from m/s and microseconds

        # compute the dt based meshgrid with z_mm displacement
        mesh_t0 = self.get_mm_coords(self.mesh['x_px'], self.mesh['y_px'], z_mm=0)
        mesh_t1 = self.get_mm_coords(self.mesh['x_px'], self.mesh['y_px'], z_mm=z_disp)

        # make sure the particles are randomly placed within the field of view of BOTH cameras
        fov = self._get_overlap_fov(mesh_t0)

        # generate a bunch of random particles within the overlaping field of view area.
        print("Creating {0} particles!".format(n_particles))
        for n in xrange(n_particles):
            p_x_mm = np.random.uniform(fov['x_mm_min'], fov['x_mm_max'])
            p_y_mm = np.random.uniform(fov['y_mm_min'], fov['y_mm_max'])
            p_z_mm = np.random.uniform(0, light_sheet_thickness) - (light_sheet_thickness / 2.0)

            self.particles_a.append(Particle(p_x_mm, p_y_mm, p_z_mm))
            self.particles_b.append(Particle(p_x_mm + x_disp, p_y_mm + y_disp, p_z_mm + z_disp))

        # now combine all of the intensity profiles for each particle into one image
        print("Projecting particles into camera image planes")
        self.images['La'], self.images['Ra'] = self._get_intensities(self.particles_a, mesh_t0, particle_size,
                                                                     particle_scatter, light_sheet_thickness)

        # and do the same thing for the displaced particles representing time = dt
        print("Displacing particles and re-projecting")
        self.images['Lb'], self.images['Rb'] = self._get_intensities(self.particles_b, mesh_t1, particle_size,
                                                                     particle_scatter, light_sheet_thickness)

        # save the images
        if output_dir is not None:
            for key in self.images.keys():
                outname = os.path.join(output_dir, "{0}{1}.tif".format(name, key))
                if os.path.exists(outname):
                    os.remove(outname)
                save_array_as_dtype(self.images[key], dtype, outname)

        # now save a json file in the same directory with the parameters of this function
        argdict = {}
        for arg in ["n_particles", "dt", "u", "v", "w", "particle_size", "particle_scatter",
                    "light_sheet_thickness", "dtype", "name", "output_dir"]:
            argdict[arg] = locals()[arg]

        with open(os.path.join(output_dir, "{0}.json".format(name)), 'w+') as logfile:
            logfile.write(json.dumps(argdict, indent=4))

        print("time = {0} hours".format((datetime.now() - start).seconds / 3600.00))


# testing area for profiling
if __name__ == '__main__':

    def main():
        mdims = (1024/4, 1280/4)
        apiv = ArtificialPIV(mdims, name="profile_test")
        apiv.load_calibration_file('cal_data/station_7/ely_may28th.cal')
        apiv.make_image_pairs(n_particles=mdims[0] * mdims[1] / 10,
                              dt=40,
                              u=5.6,
                              v=5.6,
                              w=19.0,
                              particle_size=0.2,
                              particle_scatter=100,
                              light_sheet_thickness=3.00,
                              dtype="uint16",
                              name="profile_test",
                              output_dir="artificial_images")

    import cProfile
    cProfile.run("main()")
