### PYTHON FOR MD ###
"""
    Script for merging and shifting droplet+substrate systems
    (old version: 20nm droplet)
"""

import os
from math import ceil
import mdconf as mdc

file_substrate = '/home/michele/python_for_md/FreeEnergyCorrugated/GRO/sub_height4.gro'
file_droplet = '/home/michele/python_for_md/FreeEnergyCorrugated/GRO/water_droplet.gro'
file_system = '/home/michele/python_for_md/FreeEnergyCorrugated/init_conf_h4.gro'

mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_system )

# Shifting water droplet towards the substrate
file_shifted = '/home/michele/python_for_md/FreeEnergyCorrugated/init_conf_shift_h4.gro'
mdc.shift_droplet_gro ( -4.0, 'z', file_system, file_shifted, h=1.75)

