import os
from math import ceil
import mdconf as mdc

substrates_folder = '/home/michele/python_for_md/Droplet50nmExp/Substrates'
substrates_label = [    'height1', 
                        'height2', 
                        'height3', 
                        'height4', 
                        'lambda1', 
                        'lambda2', 
                        'lambda3', 
                        'lambda4'   ]
main_folder = '/home/michele/python_for_md/Droplet50nmExp'
folders_label = [   'H1', 
                    'H2', 
                    'H3', 
                    'H4', 
                    'L1', 
                    'L2', 
                    'L3', 
                    'L4'    ]

file_droplet = '/home/michele/python_for_md/Droplet50nmExp/water_droplet.gro'

# Merge
for k in range(len(substrates_label)) :
    file_merge = main_folder+'/'+folders_label[k]+'/merge_'+substrates_label[k]+'.gro'
    file_substrate = substrates_folder+'/sub_'+substrates_label[k]+'.gro'
    mdc.merge_to_substrate_gro( file_substrate, file_droplet, file_merge )

# This will be done separately
# file_topology = '2dCylDropRoughSub/MdSystem/system_topology.top'
# mdc.add_atoms_topology(file_system, file_topology)

# Shift
delta = -17.5
height = 6.5
for k in range(len(substrates_label)) :
    file_shift = main_folder+'/'+folders_label[k]+'/shift_'+substrates_label[k]+'.gro'
    file_merge = main_folder+'/'+folders_label[k]+'/merge_'+substrates_label[k]+'.gro'
    mdc.shift_droplet_gro ( delta, 'z', file_merge, file_shift, h=height)


