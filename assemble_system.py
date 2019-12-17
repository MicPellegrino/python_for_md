import os
from math import ceil
import mdconf as mdc

file_substrate = '2dCylDropRoughSub/quad_substrate_wave_eq.pdb'
file_droplet = '2dCylDropRoughSub/droplet.pdb'
file_system = '2dCylDropRoughSub/MdSystem/initial_conf.pdb'

mdc.merge_to_substrate( file_substrate, file_droplet, file_system)

file_topology = '2dCylDropRoughSub/MdSystem/system_topology.top'
mdc.add_atoms_topology(file_system, file_topology)
