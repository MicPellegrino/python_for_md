import os
from math import ceil
import mdconf as mdc

file_substrate = '/home/michele/python_for_md/Droplet20nmExp/SubstrateCompleteMap/sub_lambda7.gro'
file_droplet = '/home/michele/python_for_md/Droplet20nmExp/water_droplet_thick.gro'
file_system = '/home/michele/python_for_md/Droplet20nmExp/SystemsThick/merged_lambda7.gro'

# mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_system )

# file_topology = '2dCylDropRoughSub/MdSystem/system_topology.top'
# mdc.add_atoms_topology(file_system, file_topology)

# Shifting water droplet towards the substrate
file_shifted = '/home/michele/python_for_md/Droplet20nmExp/SystemsThick/shifted_lambda7.gro'
mdc.shift_droplet_gro ( -5.5, 'z', file_system, file_shifted, h=1.75)

