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
# frequency = np.zeros(5, dtype=float)
# for k in range(5) :
#     frequency[k] = (0.75+0.25*k)*(1.0/amplitude)
frequency = np.zeros(12, dtype=float)
for k in range(12) :
     frequency[k] = (0.75+0.25*k)*(1.0/amplitude)

# Roughness parameter
a = np.zeros(len(frequency)+1, dtype=float)
a[0] = 0.0
for k in range(1,len(frequency)+1) :
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

# Fitting the charge-angle relation
charge_data = np.append( charges, 0.79 )
angle_data = np.append( flat_angles, 0.0 )
p = np.polyfit(charge_data, angle_data, 2)
q_val = np.linspace(charge_data[0], charge_data[-1], 250)
t_val = np.polyval(p, q_val)

# Plotting

# plt.plot(xi_fob, theta_fob, 'kx', markersize=8, label='hydrophobic')
# plt.plot(xi_fil, theta_fil, 'ks', markersize=8, label='hydrophylic')
plt.plot(xi_r, wave_angle, '.-')
plt.xlabel('Roughness parameter [nondim.]', fontsize=20.0)
plt.ylabel('Contact angle [deg]', fontsize=20.0)
# plt.legend(loc='lower left')
plt.title('Wenzel law for some reference c.a.', fontsize=20.0)
plt.show()

plt.plot(a, xi_r, 'k.-')
plt.xlabel('a', fontsize=20.0)
plt.ylabel('r', fontsize=20.0)
plt.title('Roughness param. as a function of combined height and frequency', fontsize=20.0)
plt.show()

plt.plot(q_val, t_val, 'r-')
plt.plot(q_val, t_val+3.0, 'r--')
plt.plot(q_val, t_val-3.0, 'r--')

plt.plot(charge_data, angle_data, 'ko')
plt.xlabel('Charge [e]', fontsize=20.0)
plt.ylabel('Angle [deg]', fontsize=20.0)
plt.title('Fitted charge-angle relation', fontsize=20.0)
plt.show()
