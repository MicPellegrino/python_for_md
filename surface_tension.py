"""
    SCRIPT TO GENERATE THE CONFIGURATIONS TO ESTIMATE THE SURFACE TENSION FOR 
    LIQUID-SOLID, LIQUID-VAPOUR AND SOLID-VAPOUR INTERFACES
"""

import numpy as np
import mdconf_oop as md
import os

working_folder = '/home/michele/python_for_md/SurfaceTension'
file_name_flat_SV = working_folder+'/init_SV.gro'
file_name_flat_LV = working_folder+'/init_LV.gro'
file_name_flat_SL = working_folder+'/init_SL.gro'

# Box dimensions [nm]
Lx = 10.0
Lz = Lx
Ly = Lx

# Wall position
zw = 0.5*Lz

conf_SV = md.Configuration(Lx, Ly, Lz, 'SV')
conf_LV = md.Configuration(Lx, Ly, Lz, 'LV')
conf_SL = md.Configuration(Lx, Ly, Lz, 'SL')

conf_SV.silica_monolayer(zw, cut_box=True)
conf_SV.output(file_name_flat_SV)

# Silica wall + film
file_name_flat_sol = working_folder+'/solvate.gro'
os.system('gmx solvate -cp '+file_name_flat_SV+' -cs wat_equil.gro -o '+file_name_flat_sol)
conf_SL.input(file_name_flat_sol)
conf_SL.carve_rectangle([0.0, conf_SL.box_xx], [0.25*conf_SL.box_zz, 0.75*conf_SL.box_zz])
conf_SL.output(file_name_flat_SL)

# Only water film
file_name_flat_emp = working_folder+'/empty.gro'
file_name_flat_box = working_folder+'/water_box.gro'
conf_LV.output(file_name_flat_emp)
os.system('gmx solvate -cp '+file_name_flat_emp+' -cs wat_equil.gro -o '+file_name_flat_box)
conf_LV.input(file_name_flat_box)
conf_LV.carve_rectangle([0.0, conf_LV.box_xx], [0.25*conf_LV.box_zz, 0.75*conf_LV.box_zz])
conf_LV.output(file_name_flat_LV)
