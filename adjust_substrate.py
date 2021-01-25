dx_target = 10.0            # nm
x1_upp = 0.1*1696.369873    # nm
x0_upp = 0.1*1596.369873    # nm
x1_low = -0.1*98.879997     # nm
x0_low = 0.1*1.120000       # nm

delta_x_upp = x1_upp-x0_upp
dx_upp = delta_x_upp-dx_target
print("dx_upp = ", dx_upp, "nm")
delta_x_low = x1_low-x0_low
dx_low = delta_x_low+dx_target
print("dx_low = ", dx_low, "nm")

import mdconf as md

conf_pre = 'restraints_lambda1_10nm.gro' 
conf_post = 'restraints_lambda1_10nm_adj.gro'

md.adjust_layers(dx_upp, dx_low, 15.0, conf_pre, conf_post)
