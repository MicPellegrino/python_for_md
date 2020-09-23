### PYTHON FOR MD ###
"""
    Current script for creating corrugated substartes
    (uncomment some parts!)
"""

import os
import mdconf as mdc
import numpy as np
from math import ceil
from math import sqrt
from math import pi
import scipy as sc
import scipy.special

# GROMACS version
gmx = 'gmx20'
# gmx = 'gmx'

# substrates_dir = '/home/michele/python_for_md/Droplet50nmExp/Substrates'
substrates_dir = '/home/michele/python_for_md/FreeEnergyCorrugated'

# Roughness factor
f_rough = lambda a2 : (2.0/pi) * np.sqrt(a2+1.0) * sc.special.ellipe(a2/(a2+1.0))

# Periodic images in y direction
per_y = 1.0

# Domain lenght
# Lx = 1750.0                       # [Å]
Lx = 100.0                      
Ly = per_y*46.70                    # [Å]
# Lz = 1000.0                       # [Å]
Lz = 100.0
# bf = 500.0                        # [Å]

# Lattice parameters for silica
sp = 4.50                           # [Å]
alpha_1 = sqrt(3.0/4.0)
alpha_2 = sqrt(2.0/3.0)
dy = sp*alpha_1
dz = sp*alpha_2
nj = int( np.round(Ly/(sp*alpha_1)) )
nk = 1

  ###########################################################
 ### CASE 1: fixing substrate height and change wavelenght ###
  ###########################################################

# Parameters defining the space substrate structure
h = 10.0                            # [Å]
h_off = 8*h                         # [Å]
w_off = 0.0                         # [Å^-1]
bend = True

box_y = nj*dy
box_z = nk*dz+h_off

# Flat
# ni = int( np.round((Lx+bf)/sp) )
ni = int( np.round(Lx/sp) )
file_substrate_pdb = substrates_dir+'/sub_flat.pdb'
file_uncut_gro = substrates_dir+'/uncut_flat.gro'
mdc.quad_substrate(ni, nj, nk, file_substrate_pdb)
print( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
# NB in .gro files units are expressed in nanometers
file_substrate_gro = substrates_dir+'/sub_flat.gro'
mdc.shift_and_resize_gro(-Ly/10, 0.0, 0.0, Lx/10, box_y/10, box_z/10, file_uncut_gro, file_substrate_gro)

print('\n### FLAT ###\n')
print('h = '+str(0))
print('r = '+str(1))
print('a = '+str(0))
print('omega = 0')
print('lambda = inf')

# Rough
for idx in range(4) :
    w_n = (0.5+idx*0.5)*(1.0/h)   # [Å^-1]
    a2 = (w_n*h)**2
    r0 = f_rough(a2)
    # ni = int( np.round(r0*(Lx+bf)/sp) )
    file_substrate_pdb = substrates_dir+'/sub_lambda'+str(int(idx+1))+'.pdb'
    file_uncut_gro = substrates_dir+'/uncut_lambda'+str(int(idx+1))+'.gro'
    mdc.quad_substrate_wave(h, h_off, w_n, w_off, bend, ni, nj, nk, file_substrate_pdb)
    print( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
    os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
    # NB in .gro files units are expressed in nanometers
    file_substrate_gro = substrates_dir+'/sub_lambda'+str(int(idx+1))+'.gro'
    mdc.shift_and_resize_gro(-Ly/10, 0.0, -3.0*h/10.0, Lx/10, box_y/10, box_z/10, file_uncut_gro, file_substrate_gro)
    print('\n### LAMBDA '+str(idx+1)+' ###\n')
    print('h = '+str(h))
    print('r = '+str(r0))
    print('a = '+str(np.sqrt(a2)))
    print('omega = '+str(w_n))
    print('lambda = '+str(2*np.pi/w_n))

  ###########################################################
 ### CASE 2: fixing substrate wavelenght and change height ###
  ###########################################################

# Probing a = 0.25 and a = 0.50
lambda_star = 35.0              # [Å]
w_0 = 2.0*np.pi/lambda_star     # [Å^-1]

# Flat substrate has already been produced
# Rough
for idx in range(4) :
    a = (0.5+idx*0.5)    # [nondim.]
    h_var = a/w_0
    r0 = f_rough(a**2)
    # ni = int( np.round(r0*(Lx+bf)/sp) )
    file_substrate_pdb = substrates_dir+'/sub_height'+str(int(idx+1))+'.pdb'
    file_uncut_gro = substrates_dir+'/uncut_height'+str(int(idx+1))+'.gro'
    mdc.quad_substrate_wave(h_var, h_off, w_0, w_off, bend, ni, nj, nk, file_substrate_pdb)
    print( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
    os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
    # NB in .gro files units are expressed in nanometers
    file_substrate_gro = substrates_dir+'/sub_height'+str(int(idx+1))+'.gro'
    mdc.shift_and_resize_gro(-Ly/10, 0.0, -3.0*h/10.0, Lx/10, box_y/10, box_z/10, file_uncut_gro, file_substrate_gro)

    print('\n### HEIGHT '+str(idx+1)+' ###\n')
    print('h = '+str(h_var))
    print('r = '+str(r0))
    print('a = '+str(a))
    print('omega = '+str(w_0))
    print('lambda = '+str(2*np.pi/w_n))

"""
    OLD STUFF
"""

"""
a2 = 0.50
h_var = a2/w_0
r0 = f_rough(a2)
print('\n### H | a2=0.5 ###\n')
print('h = '+str(h_var))
print('r = '+str(r0))
print('a = '+str(np.sqrt(a2)))
print('omega = '+str(w_0))
print('lambda = '+str(2*np.pi/w_0))

ni = int( np.round(r0*(Lx+bf)/sp) )
file_substrate_pdb = '/home/michele/python_for_md/Droplet20nmExp/Substrates/sub_height1.pdb'
file_uncut_gro = '/home/michele/python_for_md/Droplet20nmExp/Substrates/uncut_height1.gro'
mdc.quad_substrate_wave(h_var, h_off, w_0, w_off, bend, ni, nj, nk, file_substrate_pdb)
print( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
# NB in .gro files units are expressed in nanometers
file_substrate_gro = '/home/michele/python_for_md/Droplet20nmExp/Substrates/sub_height1.gro'
mdc.shift_and_resize_gro(-Ly/10, 0.0, 0.0, Lx/10, box_y/10, box_z/10, file_uncut_gro, file_substrate_gro)

a2 = 0.25
h_var = a2/w_0
r0 = f_rough(a2)
print('\n### H | a2=0.5 ###\n')
print('h = '+str(h_var))
print('r = '+str(r0))
print('a = '+str(np.sqrt(a2)))
print('omega = '+str(w_0))
print('lambda = '+str(2*np.pi/w_0))

ni = int( np.round(r0*(Lx+bf)/sp) )
file_substrate_pdb = '/home/michele/python_for_md/Droplet20nmExp/Substrates/sub_height2.pdb'
file_uncut_gro = '/home/michele/python_for_md/Droplet20nmExp/Substrates/uncut_height2.gro'
mdc.quad_substrate_wave(h_var, h_off, w_0, w_off, bend, ni, nj, nk, file_substrate_pdb)
print( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_uncut_gro )
# NB in .gro files units are expressed in nanometers
file_substrate_gro = '/home/michele/python_for_md/Droplet20nmExp/Substrates/sub_height2.gro'
mdc.shift_and_resize_gro(-Ly/10, 0.0, 0.0, Lx/10, box_y/10, box_z/10, file_uncut_gro, file_substrate_gro)
"""

"""
# Determine new dimensions and resize water box
file_water_in = 'wat_equil.pdb'
file_water_out = '2dCylDropRoughSubTest4/wat_equil_ext.pdb'
qs = open(file_substrate_pdb, 'r')
we = open(file_water_in, 'r')
line1_qs = qs.readline().split()
line1_we = we.readline().split()
# safety check
print(line1_qs)
print(line1_we)
alpha = 0.525
nx = int(ceil(float(line1_qs[1])/float(line1_we[1])))
# Add 1 to ensure the droplet won't shrink past lattice y dimension after imposing EOS
ny = int(ceil(float(line1_qs[2])/float(line1_we[2])))+1
nz = int(ceil(alpha*float(line1_qs[1])/float(line1_we[1])))
we.close()
qs.close()
print("gmx genconf -f " + file_water_in + " -o " + file_water_out + " -nbox %.3f %.3f %.3f" % (nx, ny, nz))
os.system("gmx genconf -f " + file_water_in + " -o " + file_water_out + " -nbox %.3f %.3f %.3f" % (nx, ny, nz))

# Carve droplet
file_droplet = '2dCylDropRoughSubTest4/droplet/wat_droplet.pdb'
beta = 0.225
mdc.carve_2D_droplet(beta, file_water_out, file_droplet, 'p')
"""
