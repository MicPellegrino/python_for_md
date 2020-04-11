import mdconf as mdc
import os
import numpy as np
from math import sqrt
from math import pi

# Finding the value of 'a'
import scipy as sc
import scipy.special
import scipy.optimize

# theta_w = 36.0
theta_w = 50.0
theta_y = 70.0

r0 = np.cos( np.deg2rad(theta_w) ) / np.cos( np.deg2rad(theta_y) )
print(r0)

f_rough = lambda a : (2.0/pi) * np.sqrt(a+1.0) * sc.special.ellipe(a/(a+1.0)) - r0
a = sc.optimize.bisect(f_rough, 4.0, 6.0)
print(a)

"""
f_rough = lambda a : (2.0/pi) * np.sqrt(a+1.0) * sc.special.ellipe(a/(a+1.0))
a = 5.0
r = f_rough(a)
print(r)
theta_y = 70.0
theta_w = np.rad2deg(np.arccos( r*np.cos(np.deg2rad(theta_y)) ))
print(theta_w)
"""
