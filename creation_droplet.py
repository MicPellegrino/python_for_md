import os
from math import ceil
import mdconf as mdc

# Determine new dimensions
# file_substrate = '2dCylDropRoughSub/quad_substrate_wave.pdb'
# file_water_in = '2dCylDropRoughSub/wat_equil.pdb'
# file_water_out = '2dCylDropRoughSub/wat_equil_ext.pdb'

file_water_in = 'wat_equil.pdb'
file_water_out = 'PreprocessingLarge/water_box.pdb'

# qs = open(file_substrate, 'r')
we = open(file_water_in, 'r')

# line1_qs = qs.readline().split()
line1_we = we.readline().split()

# alpha = 0.5
# nx = int(ceil(float(line1_qs[1])/float(line1_we[1])))
# ny = int(ceil(float(line1_qs[2])/float(line1_we[2])))
# nz = int(ceil(alpha*float(line1_qs[1])/float(line1_we[1])))

Lx = 3000.0     # [Å]
Ly = 46.70      # [Å]
Lz = 2000.0     # [Å]

nx = int( ceil( Lx/float(line1_we[1]) ) )
ny = int( ceil( Ly/float(line1_we[2]) ) )
nz = int( ceil( Lz/float(line1_we[3]) ) )

we.close()
# qs.close()

os.system("gmx20 genconf -f "+file_water_in+" -o "+file_water_out+" -nbox %.3f %.3f %.3f" % (nx, ny, nz))

file_droplet = 'PreprocessingLarge/water_droplet.pdb'

radius_factor = 0.5
mdc.carve_2D_droplet(radius_factor, file_water_out, file_droplet, 'p')
