import mdconf_oop as md

"""
conf_drop = md.Configuration(0.0, 0.0, 0.0, "Water droplet")
conf_drop.input('wat_box_3d.gro')
conf_drop.carve_sphere(7.0)
conf_drop.print_info()

conf_sub = md.Configuration(0.0, 0.0, 0.0, "Substrate")
conf_sub.input('/home/michele/python_for_md/Droplet3D_patches/substrate.gro')
conf_sub.box_zz = 23.58156
conf_sub.print_info()

conf_sub.merge(conf_drop)
conf_sub.header="Substrate + droplet system"
conf_sub.print_info()

delta_z = -0.1*(50.0-20.0)
conf_sub.shift(0.0, 0.0, delta_z)

conf_sub.output('droplet_rand_surface.gro')
"""

h = 0.3
k = 0.7/h

"""
test_conf = md.Configuration(Lx=100.0, Ly=4.67650, Lz=20.0)
test_conf.silica_monolayer(19.0, resname='WLF', cut_box=True)
test_conf.silica_monolayer_rough( 1.5, h, k, 0.0, resname='WLR' )
test_conf.output('test_oop_corrugated.gro')
"""

test_conf = md.Configuration()
test_conf.input('solv_double_layer_rough.gro')
test_conf.carve_rectangle([34.0, 66.0], [1.5+h+md.d_so, 19.0-md.d_so])
test_conf.output('test_shear_rough_init.gro')
