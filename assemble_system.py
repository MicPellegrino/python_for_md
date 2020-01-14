import os
from math import ceil
import mdconf as mdc

file_substrate = '2dCylDropRoughSubTest3/substrates/substrate_wave4_eq.gro'
file_droplet = '2dCylDropRoughSubTest3/water/droplet_eos.gro'
file_system = '2dCylDropRoughSubTest3/init_configurations/init4.gro'

mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_system )

# file_topology = '2dCylDropRoughSub/MdSystem/system_topology.top'
# mdc.add_atoms_topology(file_system, file_topology)

# Shifting water droplet towards the substrate
file_shifted = '2dCylDropRoughSubTest3/init_configurations/init4_shift.gro'
mdc.shift_droplet_gro ( -13.5, 'z', file_system, file_shifted)

# Perform equilibration
# ...

# Adding initial velocity to the water atoms in the .gro file
# file_vel = '...'
# mdc.add_velocity_gro ( 'z', -0.01, '...', file_vel )
