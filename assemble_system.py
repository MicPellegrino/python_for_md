import os
from math import ceil
import mdconf as mdc

file_substrate = '2dCylDropRoughSubTest4/substrates/equilibrated/md/wave0/confout.gro'
file_droplet = '2dCylDropRoughSubTest4/droplet/wat_droplet_eq.gro'
file_system = '2dCylDropRoughSubTest4/system_dynamic/dyn0.gro'

mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_system )

# file_topology = '2dCylDropRoughSub/MdSystem/system_topology.top'
# mdc.add_atoms_topology(file_system, file_topology)

# Shifting water droplet towards the substrate
file_shifted = '2dCylDropRoughSubTest4/system_dynamic/dyn0_shift.gro'
mdc.shift_droplet_gro ( -11.0, 'z', file_system, file_shifted)

