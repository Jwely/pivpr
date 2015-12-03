__author__ = 'Jwely'

from py.utils.get_rel_humidity import get_rel_humidity



class Experiment:

    def __init__(self, experiment_dir, experiment_id, n_samples, z_location, v_nominal, dt, test_date,
                 v_fs_mean, v_fs_sigma, q, pres_atm, temp_tunnel, wet_bulb, dry_bulb):
        """
        Class to represent the entire state of an experiment

        :param experiment_dir:  filepath to directory containing experiment folders (str)
        :param experiment_id:   integer ID number 1 through 70 (integer)
        :param n_samples:       number of observations, usually 200 (integer)
        :param z_location:      downstream distance of interrogation plane (inches, float)
        :param v_nominal:       nominal velocity target (m/s float)
        :param dt:              time between laser pulses (microseconds integer)
        :param test_date:       date on which the experiment took place (datetime)
        :param v_fs_mean:       measured free stream velocity (m/s float)
        :param v_fs_sigma:      instability in free stream velocity (m/s float)
        :param q:               stagnation pressure inside the tunnel (Pascals integer)
        :param pres_atm:        atmospheric pressure (Pascals integer)
        :param temp_tunnel:     temperature in the tunnel (Celcius float)
        :param wet_bulb:        wet bulb temperature at time of test (Kelvin float)
        :param dry_bulb:        wet bulb temperature at time of test (Kelvin float)
        :return:
        """

        # attributes as arguments (must have been recorded at time of experiment)
        self.experiment_id = experiment_id
        self.n_samples = n_samples
        self.z_location = z_location
        self.v_nominal = v_nominal
        self.dt = dt
        self.test_date = test_date
        self.v_fs_mean = v_fs_mean
        self.v_fs_sigma = v_fs_sigma
        self.q = q
        self.pres_atm = pres_atm
        self.temp_tunnel = temp_tunnel
        self.wet_bulb = wet_bulb
        self.dry_bulb = dry_bulb

        # calculate relative humidity from inputs
        self.rel_humid = get_rel_humidity(dry_bulb, wet_bulb, pres_atm)

        # attributes ingested by function call
        self.axial_vortex = None        # will be an actual AxialVortex instance

        # calculated and derived attributes from the axial vortex
        self.core_radius = None         # distance between core and vmax_tangential location (mm)
        self.vmax_tangential = None     # maximum tangential velocity
        self.vmax_axial = None          # maximum axial velocity
        self.vcore_axial = None         # axial velocity at the core


    def ingest_axial_vortex(self, axial_vortex_instance):
        """
        :param axial_vortex_instance:
        :return:
        """
        self.axial_vortex = axial_vortex_instance