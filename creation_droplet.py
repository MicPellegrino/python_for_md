## PYTHON FOR MD ##
"""
    Script for carving water droplets out of a single block of SCP/E 
    pre-equilibrated water
"""

import os
from math import ceil
import mdconf as mdc

# Determine new dimensions
# file_substrate = '2dCylDropRoughSub/quad_substrate_wave.pdb'
# file_water_in = '2dCylDropRoughSub/wat_equil.pdb'
# file_water_out = '2dCylDropRoughSub/wat_equil_ext.pdb'

file_water_in = 'wat_equil.pdb'
file_water_out = '/home/michele/python_for_md/FreeEnergyCorrugated/water_box_thick.pdb'

# qs = open(file_substrate, 'r')
we = open(file_water_in, 'r')

# line1_qs = qs.readline().split()
line1_we = we.readline().split()

# alpha = 0.5
# nx = int(ceil(float(line1_qs[1])/float(line1_we[1])))
# ny = int(ceil(float(line1_qs[2])/float(line1_we[2])))
# nz = int(ceil(alpha*float(line1_qs[1])/float(line1_we[1])))

# Large droplet configuration
"""
Lx = 3000.0     # [Å]
Ly = 46.70      # [Å]
Lz = 2000.0     # [Å]
"""
# Medium droplet configuration; longer periodic dimension
Lx = 750.000        # [Å]
Ly = 23.3830        # [Å]
Lz = 353.7240       # [Å]

nx = int( ceil( Lx/float(line1_we[1]) ) )
ny = int( ceil( Ly/float(line1_we[2]) ) )
nz = int( ceil( Lz/float(line1_we[3]) ) )

we.close()
# qs.close()

os.system("gmx20 genconf -f "+file_water_in+" -o "+file_water_out+" -nbox %.3f %.3f %.3f" % (nx, ny, nz))

file_droplet = '/home/michele/python_for_md/FreeEnergyCorrugated/water_droplet.pdb'

radius = 100.0      #[Å]
radius_factor = 2.0*radius/Lz

mdc.carve_2D_droplet(radius_factor, file_water_out, file_droplet, 'p')
