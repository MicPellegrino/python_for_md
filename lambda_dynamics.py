import numpy as np

# Substrate difference PER LAYER (i.e. upper or lower)
DL_x = 0.5              # [nm]

# Desired shear rate
# gamma = 1.0e-3/48.0   # [ps^-1]

# Distance between bilayers
<<<<<<< HEAD
<<<<<<< HEAD
# L_z = 29.0919983        # [nm]
L_z = 9.933             # [nm]
=======
# L_z = 29.0919983          # [nm]
L_z = 2.3299999-0.385000    # [nm]
>>>>>>> 569e7de7324a4fa570dd3140ab44c5d11bf41660
=======
L_z = 29.0919983          # [nm]
# L_z = 2.3299999-0.385000    # [nm]
# L_z = 50.0
>>>>>>> d539aed3ad013b3023661f1d43045f8e4ea15547

# Total velocity (upper-lower)
# v_tot = gamma*L_z     # [nm/ps]

# Velocity PER LAYER (i.e upper or lowe, wrt initial restraints)
# v_lay = 0.5*v_tot     # [nm/ps]
<<<<<<< HEAD
v_lay = 10.0e-3      # [nm/ps]
=======
v_lay = 3.33e-3         # [nm/ps]
>>>>>>> 569e7de7324a4fa570dd3140ab44c5d11bf41660

# Time step of the simulation
dt = 0.004              # [ps]

# Time difference between lambda=1 and lambda=0
Dtime = DL_x/v_lay      # [ps]

# Number of steps to increase lambda
N_lambda = Dtime/dt     # [nondim.]
print("Nl = "+str(N_lambda))

# Increment in lambda at each step
d_lambda = 1.0/N_lambda
print("dl = "+str(d_lambda))
