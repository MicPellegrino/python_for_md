import os
from math import ceil
import mdconf as mdc

# Testing
# os.system("gmx help genconf")

# Create flat substrate
file_flat_substrate = '2dCylDropRoughSubTest3/quad_substrate_flat.pdb'
ni = 150
nj = 20
nk = 1
# mdc.quad_substrate(ni, nj, nk, file_flat_substrate)

# Only one substrate, for now

h = 4.0
h_off = 2*h
w_off = 0.0
bend = True

for idx in range(5) :
    w_n = (0.75+idx*0.25)*(1.0/h)
    file_substrate = '2dCylDropRoughSubTest3/quad_substrate_wave'+str(idx)+'.pdb'
    # mdc.quad_substrate_wave(h, h_off, w_n, w_off, bend, ni, nj, nk, file_substrate)

# Determine new dimensions and resize wateer box
file_water_in = '2dCylDropRoughSubTest3/wat_equil.pdb'
file_water_out = '2dCylDropRoughSubTest3/wat_equil_ext.pdb'
qs = open(file_substrate, 'r')
we = open(file_water_in, 'r')
line1_qs = qs.readline().split()
line1_we = we.readline().split()
print(line1_qs)
print(line1_we)
alpha = 0.525
nx = int(ceil(float(line1_qs[1])/float(line1_we[1])))
# Add 1 to ensure the droplet won't shrink past lattice y dimension after imposing EOS
ny = int(ceil(float(line1_qs[2])/float(line1_we[2])))+1
nz = int(ceil(alpha*float(line1_qs[1])/float(line1_we[1])))
we.close()
qs.close()
print("gmx genconf -f "+file_water_in+" -o "+file_water_out+" -nbox %.3f %.3f %.3f" % (nx, ny, nz))
os.system("gmx genconf -f "+file_water_in+" -o "+file_water_out+" -nbox %.3f %.3f %.3f" % (nx, ny, nz))

# Carve droplet
file_droplet = '2dCylDropRoughSubTest3/wat_droplet.pdb'
beta = 0.125
mdc.carve_2D_droplet(beta, file_water_out, file_droplet, 'p')

# Equilibrate droplet
# ...

# EOS droplet
# ...

# Equilibrate substrates
# NB! RESIZE SYSTEM HEIGHT SO THAT EWALD SUMMATION CAN BE PERFORMED!
# ...

# file_system = '2dCylDropRoughSubTest3/system.pdb'
# mdc.merge_to_substrate( file_substrate, file_droplet, file_system)
