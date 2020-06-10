import numpy as np

# Substrate difference PER LAYER (i.e. upper or lower)
DL_x = 1.0              # [nm]

# Desired shear rate
# gamma = 1.0e-3/48.0   # [ps^-1]

# Distance between bilayers
L_z = 29.0919983        # [nm]

# Total velocity (upper-lower)
# v_tot = gamma*L_z     # [nm/ps]

# Velocity PER LAYER (i.e upper or lowe, wrt initial restraints)
# v_lay = 0.5*v_tot     # [nm/ps]
v_lay = 10.0e-3      # [nm/ps]

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
