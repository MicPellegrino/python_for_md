import mdconf_oop as md
import os
import numpy as np

"""
    Scrip to produce several pos-res configurations of corrugate substrate
"""

working_folder = '/home/michele/python_for_md/FreeEnergyCorrugate'
file_name_flat_res = working_folder+'/flat_res.gro'
file_name_flat_sol = working_folder+'/flat_sol.gro'
file_name_flat_init = working_folder+'/flat_init.gro'

# Nominal lenght along x
Lx_nominal = 20.0
# Nominal lenght along z
Lz_nominal = 10.0

conf_flat = md.Configuration(Lx=Lx_nominal, Ly=2.33827, Lz=Lz_nominal)

zw = 5.0
conf_flat.silica_monolayer(zw, cut_box=True)
conf_flat.output(file_name_flat_res)

n_quad_x = int( conf_flat.box_xx/0.450 )

# Fix wave number
k = 2*np.pi/(conf_flat.box_xx/4)
a_vec = np.linspace(0.1, 1.5, 15)
h_vec = a_vec/k
labels = np.linspace(0, len(h_vec)-1, len(h_vec)).astype(int)
conf_rough_vec = []

area_const = conf_flat.box_xx * conf_flat.box_zz

for l in range(len(labels)) :
    conf_rough_vec.append(md.Configuration(Lx=Lx_nominal, Ly=2.33827, Lz=Lz_nominal))
    conf_rough_vec[l].silica_monolayer_rough(zw, h_vec[l], k, 0.5*np.pi, mode='number', ni=n_quad_x, cut_box=True)
    Lzk = area_const / conf_rough_vec[l].box_xx
    conf_rough_vec[l].box_zz = Lzk
    conf_rough_vec[l].shift(0.0, 0.0, 0.5*(Lzk-Lz_nominal), 'SUB')
    conf_rough_vec[l].output(working_folder+'/a'+str(labels[l]).zfill(2)+'_res.gro')

# Solvate
os.system('gmx solvate -cp '+file_name_flat_res+' -cs wat_equil.gro -scale 0.6 -o '+file_name_flat_sol)

# Resize
conf_init = md.Configuration()
conf_init.input(file_name_flat_sol)
conf_init.carve_rectangle([0.0, conf_flat.box_xx], [2.5, 7.5])
conf_init.output(file_name_flat_init)

conf_init.print_info()
