from math import pi
import scipy as sc
import scipy.special
import numpy as np
import matplotlib.pyplot as plt

# a = amplitude^2 * frequency^2
rough_parameter = lambda a : (2.0/pi) * np.sqrt(a+1.0) * sc.special.ellipe(a/(a+1.0))

"""
q_O (e)     theta_0 (deg.) (approx)
-0.40       110
-0.60       90
-0.67       70
-0.74       37
-0.79       0
The case on complete wetting is not considered
"""

# Substrate partial charges and respective contact angles
charges = np.array([0.40, 0.60, 0.67, 0.74])
flat_angles = np.array([110, 90, 70, 37])

# Substrate amplitude and wavelenght
amplitude = 4.0
frequency = np.zeros(5, dtype=float)
for k in range(5) :
    frequency[k] = (0.75+0.25*k)*(1.0/amplitude)

# Roughness parameter
a = np.zeros(6, dtype=float)
a[0] = 0.0
for k in range(1,6) :
    a[k] = (0.75+0.25*k)**2
xi_r = rough_parameter(a)

# Storage for contact angles
cos_wave_angle = np.outer( xi_r, np.cos( np.deg2rad(flat_angles) ) )
sat_arccos = lambda cval : 0.0 if (cval>1.0) else np.rad2deg(np.arccos(cval))
wave_angle = np.vectorize(sat_arccos)(cos_wave_angle)

# Some data
a_fob = np.array( [1.5625, 3.0625] )
xi_fob = rough_parameter(a_fob)
theta_fob = [122.01, 114.27]
a_fil = np.array( [1.5625, 3.0625] )
xi_fil = rough_parameter(a_fil)
theta_fil = [81.66, 66.35]

# Plotting
plt.plot(xi_fob, theta_fob, 'kx', markersize=8, label='hydrophobic')
plt.plot(xi_fil, theta_fil, 'ks', markersize=8, label='hydrophylic')
plt.plot(xi_r, wave_angle, 'o--')
plt.xlabel('Roughness parameter [nondim.]')
plt.ylabel('Contact angle [deg]')
plt.legend(loc='lower left')
plt.show()

