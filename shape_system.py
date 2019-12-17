import os
from math import ceil
import mdconf as mdc

# Testing
os.system("gmx help genconf")

# Determine new dimensions

file_substrate = 'md_system/quad_substrate_wave.pdb'
file_water_in = 'md_system/wat_equil_eos.pdb'
file_water_out = 'md_system/wat_equil_eos_ext.pdb'

qs = open(file_substrate, 'r')
we = open(file_water_in, 'r')

line1_qs = qs.readline().split()
line1_we = we.readline().split()

print(line1_qs)
print(line1_we)

alpha = 0.5
nx = int(ceil(float(line1_qs[1])/float(line1_we[1])))
ny = int(ceil(float(line1_qs[2])/float(line1_we[2])))
nz = int(ceil(alpha*float(line1_qs[1])/float(line1_we[1])))

we.close()
qs.close()

print("gmx genconf -f "+file_water_in+" -o "+file_water_out+" -nbox %.3f %.3f %.3f" % (nx, ny, nz))
os.system("gmx genconf -f "+file_water_in+" -o "+file_water_out+" -nbox %.3f %.3f %.3f" % (nx, ny, nz))

file_droplet = 'md_system/wat_droplet.pdb'

mdc.carve_2D_droplet(0.50, file_water_out, file_droplet, 'p')

file_system = 'md_system/system.pdb'

mdc.merge_to_substrate( file_substrate, file_droplet, file_system)
