### PYTHON FOR MD ###
"""
    OLD script for creating substrates
"""

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
theta_w = 55.0
theta_y = 70.0

r0 = np.cos( np.deg2rad(theta_w) ) / np.cos( np.deg2rad(theta_y) )
f_rough = lambda a : (2.0/pi) * np.sqrt(a+1.0) * sc.special.ellipe(a/(a+1.0)) - r0
a = sc.optimize.bisect(f_rough, 3.5, 6.5)

d_hb = 3.0                                  # [Å]
c = sqrt(a)/(2.0*pi)                        # Nondim.
# If PI_1 gets to small, the curvature of the wavy substrate may beak LINCS algorithm
# Keep it large enough, but still around one o.o.m. larger than hydrogen bond distance
# (so that interesting things may be observed...)
PI_1 = 25                                   # Nondim.
l = PI_1*d_hb                               # [Å]
h = c * l                                   # [Å]
omega = 2*pi/l                              # [Å^-1]
d_R0 = 1000                                 # [Å]
PI_2 = l/d_R0
print("PI_2 = "+str(PI_2))

Lx = 3000.0                                 # [Å]
Ly = 46.70                                  # [Å]
Lz = 2000.0                                 # [Å]

buffer = 250.0                              # [Å]

sp = 4.50                                   # [Å]
alpha_1 = sqrt(3.0/4.0)
alpha_2 = sqrt(2.0/3.0)
dy = sp*alpha_1
dz = sp*alpha_2

print('r0 = \t\t%f' % r0)
print('a = \t\t%f' % a)
print('l = \t\t'+str(l))
print('h = \t\t'+str(h))
print('omega = \t'+str(omega))

# TILTING W.R.T. WETTING FRONT
theta = 25.0

ni_f = int(np.round( Lx/sp ))
ni_w = int(np.round( r0*(Lx+buffer)/sp ))
nj = int(np.round( Ly/(sp*alpha_1*np.cos(np.deg2rad(theta))) ))
nk = 1

box_y = nj*dy
box_z = nk*dz+2*h+3.0

# CREATE WAVY SUBSTRATE (NB! .pbd format is in Å !!!)
gmx = 'gmx'
# gmx = 'gmx20'
work_dir = os.getcwd()+'/PreprocessingLarge/substrate_tilted'
fn1 = 'quad_sub_flat'
fn2 = 'quad_sub_wave'
bend = True
mdc.quad_substrate(ni_f, nj, nk, work_dir+'/'+fn1+'.pdb')
mdc.quad_substrate_wave_tilted(h, h, omega, 0.0, bend, theta, ni_w, nj, nk, work_dir+'/'+fn2+'.pdb')
os.system(gmx+' editconf -f '+work_dir+'/'+fn1+'.pdb -o '+work_dir+'/'+fn1+'.gro')
os.system(gmx+' editconf -f '+work_dir+'/'+fn2+'.pdb -o '+work_dir+'/'+fn2+'.gro')
# mdc.shift_and_resize_gro(-Ly, 0.0, 0.0, 300.0, box_y, box_z, work_dir+'/'+fn1+'.gro', work_dir+'/'+fn1+'_res.gro')
mdc.shift_and_resize_gro(-Ly/10, 0.0, 0.0, Lx/10, box_y/10, box_z/10, work_dir+'/'+fn2+'.gro', work_dir+'/'+fn2+'_res.gro')

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
