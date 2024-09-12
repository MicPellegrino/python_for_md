import numpy as np
import matplotlib.pyplot as plt
import mdconf_oop as md

def read_xvg_file(file_name) :
    x = []
    y = []
    with open(file_name) as f:
        for line in f:
            cols = line.split()
            if len(cols) == 2 and cols[0][0] != '@':
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return x, y

# Water-vapourt surface tension
gamma = 5.78e-2         # Pa*m
# Avogadro's number
N_A = 6.02214076e23     # mol^-1

# Box lenghts [m]
Lx0 = 20.70000e-9
Lx1 = 18.59721e-9
Ly = 4.67654e-9

theta_0 = 37.8
delta_F_wet_dry = -2.0*N_A*gamma*np.cos(np.deg2rad(theta_0))*Ly*Lx0*1e-3
delta_F_surface =  2.0*N_A*gamma*Ly*(Lx1-Lx0)*1e-3

xvg_file_dry = 'FreeEnergyCorrugate/FE_connect/dhdl.xvg'
xvg_file_wet = 'FreeEnergyCorrugate/FE_consolv/dhdl.xvg'

t, _ = read_xvg_file( xvg_file_dry )
lam = np.array(t) / max(t)
dl = lam[1]-lam[0]

_, dhdl_avg_dry = read_xvg_file( xvg_file_dry )
delta_h_dry = dl*np.cumsum(dhdl_avg_dry)

_, dhdl_avg_wet = read_xvg_file( xvg_file_wet )
delta_h_wet = dl*np.cumsum(dhdl_avg_wet)

# plt.plot(t, dhdl_avg_wet, 'b-', linewidth=0.25, label='wet')
# plt.plot(t, dhdl_avg_dry, 'r-', linewidth=0.25, label='dry')
plt.plot(t, delta_h_dry, 'k--', linewidth=2.0)
plt.plot(t, delta_h_wet, 'k--', linewidth=2.0)
plt.show()

print("### Free-energy integration results ###")
print("dA_dry  = "+str(delta_h_dry[-1])+" kJ/mol")
print("dA_wet  = "+str(delta_h_wet[-1])+" kJ/mol")
print("dA_tens = "+str(delta_F_surface)+" kJ/mol")
print("dA_exc  = "+str(delta_h_wet[-1]-delta_h_dry[-1]-delta_F_surface)+" kJ/mol")
print("dA_ref  = "+str(delta_F_wet_dry)+" kJ/mol")
print("#######################################")
