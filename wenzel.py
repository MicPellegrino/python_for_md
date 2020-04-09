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

# Cassie
sigma = np.linspace(0.0,1.0,100)
cos_cassie_angle = np.outer( sigma, np.cos( np.deg2rad(flat_angles) ) ) - \
    np.outer( (1.0-sigma), np.ones(len(flat_angles)) )
cassie_angle = np.vectorize(sat_arccos)(cos_cassie_angle)

# Plotting

# plt.plot(xi_fob, theta_fob, 'kx', markersize=8, label='hydrophobic')
# plt.plot(xi_fil, theta_fil, 'ks', markersize=8, label='hydrophylic')
# [p1, p2, p3, p4] = plt.plot(xi_r, wave_angle, '-')
# plt.xlabel('Roughness parameter [nondim.]', fontsize=20.0)
# plt.ylabel('Contact angle [deg]', fontsize=20.0)
# plt.legend([p1, p2, p3, p4], [r'$\theta_Y=110$', r'$\theta_Y=90$', r'$\theta_Y=70$', r'$\theta_Y=37$'], loc='upper left', fontsize=20.0)
# plt.title('Wenzel law for some reference contact angle', fontsize=20.0)
# plt.xlim([min(xi_r), max(xi_r)])
# plt.xticks(fontsize=20.0)
# plt.yticks(fontsize=20.0)
# plt.show()

[p1, p2, p3, p4] = plt.plot(sigma, cassie_angle, '-')
plt.xlabel('Solid surface fraction [nondim.]', fontsize=20.0)
plt.ylabel('Contact angle [deg]', fontsize=20.0)
plt.legend([p1, p2, p3, p4], [r'$\theta_Y=110$', r'$\theta_Y=90$', r'$\theta_Y=70$', r'$\theta_Y=37$'], loc='lower left', fontsize=20.0)
plt.title('Cassie-Baxter law for some reference contact angle', fontsize=20.0)
plt.xlim([min(sigma), max(sigma)])
plt.xticks(fontsize=20.0)
plt.yticks(fontsize=20.0)
plt.show()

# plt.plot(a, xi_r, 'k-')
# plt.xlabel('a [nondim.]', fontsize=20.0)
# plt.ylabel('r [nondim.]', fontsize=20.0)
# plt.title('Roughness paramameter as a function of combined height and frequency', fontsize=20.0)
# plt.xlim([min(a), max(a)])
# plt.ylim([min(xi_r), max(xi_r)])
# plt.xticks(fontsize=20.0)
# plt.yticks(fontsize=20.0)
# plt.show()

# plt.plot(q_val, t_val, 'r-')
# plt.plot(q_val, t_val+3.0, 'r--')
# plt.plot(q_val, t_val-3.0, 'r--')
# plt.plot(charge_data, angle_data, 'ko')
# plt.xlabel('Charge [e]', fontsize=20.0)
# plt.ylabel('Angle [deg]', fontsize=20.0)
# plt.title('Fitted charge-angle relation', fontsize=20.0)
# plt.xticks(fontsize=20.0)
# plt.yticks(fontsize=20.0)
# plt.show()
