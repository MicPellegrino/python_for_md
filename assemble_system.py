import os
from math import ceil
import mdconf as mdc

file_substrate = '2dCylDropRoughSubTest3/substrates/substrate_wave4_eq.pdb'
file_droplet = '2dCylDropRoughSubTest3/water/droplet_eos.pdb'
file_system = '2dCylDropRoughSubTest3/init_configurations/init4.pdb'

mdc.merge_to_substrate( file_substrate, file_droplet, file_system )

# file_topology = '2dCylDropRoughSub/MdSystem/system_topology.top'
# mdc.add_atoms_topology(file_system, file_topology)

# Shifting water droplet towards the substrate
file_shifted = '...'
mdc.shift_droplet ( -135, 'z', file_system, file_shifted)

# Perform equilibration
# ...

# Adding initial velocity to the water atoms in the .gro file
file_vel = '...'
mdc.add_velocity_gro ( 'z', -0.01, '...', file_vel )
