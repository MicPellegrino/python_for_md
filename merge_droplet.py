import os
from math import ceil
import mdconf as mdc

file_droplet = '/home/michele/python_for_md/Droplet20nmExp/water_droplet.gro'

# Flat
file_substrate = '/home/michele/python_for_md/Droplet20nmExp/Substrates/sub_flat.gro'
file_system = '/home/michele/python_for_md/Droplet20nmExp/Systems/system_flat.gro'
mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_system )

# Rough
for idx in range(5) :
    file_substrate = '/home/michele/python_for_md/Droplet20nmExp/Substrates/sub_wave'+str(int(idx+1))+'.gro'
    file_system = '/home/michele/python_for_md/Droplet20nmExp/Systems/system_wave'+str(int(idx+1))+'.gro'
    mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_system )
