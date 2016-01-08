__author__ = 'Jwely'




class ArtificialParticles:
    def __init__(self):


        self.dt = None      # time between snapshots

        # dictionary of actual calibration information
        self.cal_info = {"l_x_mm": [],
                         "l_y_mm": [],
                         "l_x_px": [],
                         "l_Y_px": [],
                         "r_x_mm": [],
                         "r_y_mm": [],
                         "r_x_px": [],
                         "r_y_px": []}

        # simple dictionary to translate text in the cal file into tags
        self.tag_dict = {"Left Camera x mm": "l_x_mm",
                         "Left Camera y mm": "l_y_mm",
                         "Left Camera X pix": "l_x_px",
                         "Left Camera Y pix": "l_Y_px",
                         "Right Camera x mm": "r_x_mm",
                         "Right Camera y mm": "r_y_mm",
                         "Right Camera X pix": "r_x_px",
                         "Right Camera Y pix": "r_y_px"}

    def load_calibration_file(self, filepath):
        """
        Parses a .cal file of format 'TSI_CAL_VERSION 2.100000' and extracts the equation
        information. This is expressed as a list of 9 coefficients and the corresponding orders of the
        x,y,z displacement terms to be multiplied. So a row of 123, 0, 1, 1 would translate to '123 * y * z'

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

                if line[0] in ["L", "R", "C"]:
                    active = get_active(line)
                    continue

                if active is not None:
                    coeff, xord, yord, zord = line.replace('\n', '').split(', ')
                    self.cal_info[active].append([float(coeff), int(xord), int(yord), int(zord)])

        print self.cal_info



if __name__ == "__main__":

    ap = ArtificialParticles()
    ap.load_calibration_file(r"D:\pivpr_raw_data\Calibration\Location 3\ely_may28th.cal")




