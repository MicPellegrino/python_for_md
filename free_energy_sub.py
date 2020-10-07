### PYTHON FOR MD ###
"""
    Creates the substrate configurations for free-energy calculations
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
gmx = 'gmx'

substrates_dir = '/home/michele/python_for_md/FreeEnergyCorrugated'

# Roughness factor
f_rough = lambda a2 : (2.0/pi) * np.sqrt(a2+1.0) * sc.special.ellipe(a2/(a2+1.0))

# Domain lenght
Lx = 250.0              # [Å]
Ly = 46.70              # [Å]
Lz = 100.0              # [Å]

# Lattice parameters for silica
sp = 4.50               # [Å]
alpha_1 = sqrt(3.0/4.0)
alpha_2 = sqrt(2.0/3.0)
dy = sp*alpha_1
dz = sp*alpha_2
ni = int( np.round(Lx/sp) )
nj = int( np.round(Ly/(sp*alpha_1)) )
nk = 1

h = 10.0                            # [Å]
h_off = 3*h                         # [Å]
w_off = 0.0                         # [Å^-1]

file_substrate_pdb = substrates_dir+'/sub_flat.pdb'
file_substrate_gro = substrates_dir+'/sub_flat.gro'
mdc.quad_substrate_wave(0.0, h_off, 0.0, w_off, False, ni, nj, nk, file_substrate_pdb)
os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_substrate_gro)

bend = True

# Rough
for idx in range(4) :
    w_n = (0.5+idx*0.5)*(1.0/h)   # [Å^-1]
    a2 = (w_n*h)**2
    r0 = f_rough(a2)
    file_substrate_pdb = substrates_dir+'/sub_lambda'+str(int(idx+1))+'.pdb'
    file_substrate_gro = substrates_dir+'/sub_lambda'+str(int(idx+1))+'.gro'
    mdc.quad_substrate_wave(h, h_off, w_n, w_off, bend, ni, nj, nk, file_substrate_pdb)
    os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_substrate_gro )

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
    file_substrate_gro = substrates_dir+'/sub_height'+str(int(idx+1))+'.gro'
    mdc.quad_substrate_wave(h_var, h_off, w_0, w_off, bend, ni, nj, nk, file_substrate_pdb)
    os.system( gmx+" editconf -f " + file_substrate_pdb + " -o " + file_substrate_gro )
