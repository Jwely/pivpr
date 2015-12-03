__author__ = 'Jwely'

import math


def get_rel_humidity(dry_bulb, wet_bulb, pressure):
    """
    calculates relative humidity as a function of inputs. Equations taken from
    [https://www.easycalculation.com/weather/dewpoint-humidity-calculator.php]

    :param wet_bulb:    temperature of wet bulb in K
    :param dry_bulb:    temperature of dry bulb in K
    :param pressure:    air pressure in Pascals
    :return:
    """

    t_dry = (float(dry_bulb) - 273.15)  # convert to C
    t_wet = (float(wet_bulb) - 273.15)  # convert to C
    p = float(pressure) / 100           # convert millibars

    es = 6.112 * math.exp(17.67 * t_dry / (t_dry + 243.5))
    ew = 6.112 * math.exp(17.67 * t_wet / (t_wet + 243.5))
    e = ew - p * (t_dry - t_wet) * 0.00066 * (1 + (0.00115 * t_wet))
    rh = 100 * e / es
    return rh


if __name__ == "__main__":

    print get_rel_humidity(299.15, 293.65, 102000)