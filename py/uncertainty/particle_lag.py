from math import exp

"""
scratch script to estimate particle lag
https://engineering.purdue.edu/~yanchen/paper/2014-1.pdf
"""

# params and constants
d_p = 6e-7      # particle diameter (m)
dp2 = d_p ** 2  # particle diameter squared (m ^ 2)
rho_p = 838     # density of particle (kg / m^3)
rho_f = 1.225   # density of air (kg / m^3)
mu = 1.983e-5   # dynamic viscosity of air (kg / m * s) or (pascal * s)

u = 1          # characteristic flow velocity m/s (maximum)

char_dim = rho_p * dp2 / (18 * mu) * (u / 0.1) * 1000
print(str(char_dim) + "mm")

