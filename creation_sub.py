import mdconf as mdc
import os
import numpy as np
from math import sqrt
from math import pi

# Finding the value of 'a'
import scipy as sc
import scipy.special
import scipy.optimize
theta_w = 36.0
theta_y = 70.0
r0 = np.cos( np.deg2rad(theta_w) ) / np.cos( np.deg2rad(theta_y) )
f_rough = lambda a : (2.0/pi) * np.sqrt(a+1.0) * sc.special.ellipe(a/(a+1.0)) - r0
a = sc.optimize.bisect(f_rough, 9.25, 11.85)

d_hb = 0.3                                  # [nm]
c = sqrt(a)/(2.0*pi)                        # Nondim.
l = d_hb*np.array([1.000, 10.00, 100.0])    # [nm]
h = c * l                                   # [nm]
omega = 2*pi/l                              # [nm^-1]

Lx = 150.0                                  # [nm]
Ly = 4.670                                  # [nm]

sp = 0.450                              #    [nm]
alpha_1 = sqrt(3.0/4.0)

print('r0 = \t\t%f' % r0)
print('a = \t\t%f' % a)
print('l = \t\t'+str(l))
print('h = \t\t'+str(h))
print('omega = \t'+str(omega))

# Those have to be revised!
# ni = int(np.round(Lx/sp))
# nj = int(np.round(Ly/(sp*alpha_1)))
# nk = 1

# CREATE WAVY SUBSTRATE (NB! .pbd format is in Å !!!)
# work_dir = os.getcwd()+'/substrate_large'
# fn1 = 'quad_sub_001.pdb'
# fn2 = 'quad_sub_010.pdb'
# fn3 = 'quad_sub_100.pdb'
# bend = True
# mdc.quad_substrate_wave(10*h[0], 20.0, 0.1*omega[0], 0.0, bend, ni, nj, nk, work_dir+'/'+fn1)
# mdc.quad_substrate_wave(10*h[1], 20.0, 0.1*omega[1], 0.0, bend, ni, nj, nk, work_dir+'/'+fn2)
# mdc.quad_substrate_wave(10*h[2], 20.0, 0.1*omega[2], 0.0, bend, ni, nj, nk, work_dir+'/'+fn3)

###############################################################################
### GENERATES FLAT LJ SUBSTRATE ###############################################
###############################################################################

# ni = 32
# nj = 17
# nk = 3
# sp = 2.7
# alpha_1 = sqrt(3.0/4.0)
# alpha_2 = sqrt(2.0/3.0)
# file_name = '../layering/lj_substrate.pdb'
# mdc.lj_substrate ( ni, nj, nk, sp, alpha_1, alpha_2, file_name )

###############################################################################
### GENERATES FLAT SIO2 SUBSTRATE #############################################
###############################################################################

# ni = 100
# nj = 20
# nk = 1
# file_name = '/home/michele/python_for_md/2dCylDropFlatHyst1/quad_substrate.pdb'
# mdc.quad_substrate( ni, nj, nk, file_name )

###############################################################################
### GENERATES WAVY SIO2 SUBSTRATE #############################################
###############################################################################

# h = 4.5
# h_off = h
# w_n = 1.0/h
# w_off = 0.0
# ni = 100
# nj = 20
# nk = 1
# file_name = '2dCylDropRoughSub/quad_substrate_wave.pdb'
# bend = True

# mdc.quad_substrate_wave(h, h_off, w_n, w_off, bend, ni, nj, nk, file_name)

# mdc.add_atoms_topology('energy_minimization/quad_substrate_wave.pdb', 'energy_minimization/quad_substrate_wave.top')
# mdc.add_atoms_topology('md_system/system.pdb', 'md_system/system.top')
# mdc.add_atoms_topology('md_system/equil_droplet/wat_droplet.pdb', 'md_system/equil_droplet/wat_droplet.top')
# mdc.add_atoms_topology('md_system/md_equil/system_eq.pdb', 'md_system/md_equil/system_eq.top')
