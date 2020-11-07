import numpy as np

L = 29.092e-9       # [m]
U = 0.01e3          # [m/s]
l_s = 1.0e-9        # [m]
mu = 0.0008900      # [Pa*s]
rho = 997.05        # [kg/m^3]

a = (1.0/rho) * (2.0*mu*U) / (0.25*L**2 + l_s*L)    # [m/s^2]
a = a*1e9/(1e24)

print("Acceleration:")
print("a = "+str(a)+" nm/ps^2")
