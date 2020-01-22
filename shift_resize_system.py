import mdconf as mdc

bx = 47.72970
by = 5.51140
bz = 35.37240
new_box_x = 0.5*bx
new_box_y = by
new_box_z = 0.5*bz
delta_x = -0.25*bx
delta_y = 0.0
delta_z = 0.0
input_file = 'resize_in.gro'

mdc.shift_and_resize_gro(delta_x, delta_y, delta_z, new_box_x, new_box_y, new_box_z, input_file)

