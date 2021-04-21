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

xvg_res =      [ '/home/michele/python_for_md/FreeEnergyCorrugated/FE06_res/dhdl.xvg',      ]
xvg_sol =      [ '/home/michele/python_for_md/FreeEnergyCorrugated/FE06_sol/dhdl.xvg',      ]
gro_init_res = [ '/home/michele/python_for_md/FreeEnergyCorrugated/FE06_res/init_conf.gro', ]
gro_fin_res =  [ '/home/michele/python_for_md/FreeEnergyCorrugated/FE06_res/confout.gro',   ]

conf_init = md.Configuration()
conf_init.input(gro_init_res[0])
conf_fin = md.Configuration()
conf_fin.input(gro_fin_res[0])

Ly = (conf_init.box_yy)*1e-9
Lx0 = (conf_init.box_xx)*1e-9
Lx1 = (conf_fin.box_xx)*1e-9

# Free-energy difference due to surface tension
diff_H_sur = (2.0*N_A*gamma*Ly*(Lx1-Lx0))*1e-3
print("diff_H_sur = "+str(diff_H_sur)+" kJ/mol")

t, _ = read_xvg_file( xvg_res[0] )

lam = np.array(t) / max(t)
dt = lam[1]-lam[0]
delta_lam = 2.0e-6

dhdl = np.zeros( (len(xvg_res), len(t)), dtype=float )
for k in range(len(xvg_res)) :
    _, y = read_xvg_file( xvg_res[k] )
    dhdl[k,:] = np.array(y)
dhdl_avg = np.mean(dhdl, axis=0)
delta_h = dt*np.cumsum(dhdl_avg)
diff_H_res = np.mean(dhdl_avg)
print("diff_H_res = "+str(diff_H_res)+" kJ/mol")

fig1, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(lam, dhdl_avg, 'g-', linewidth=0.25)
ax1.plot(lam, delta_h, 'r-', linewidth=2.5, label='integral')
ax1.plot(lam, np.zeros(len(lam)), 'k--', linewidth=2.5)

ax1.set_title('Dry')
ax1.set_ylabel(r'$<dH/d\lambda>$')
ax1.set_xlabel(r'$\lambda$')
ax1.legend()

max_y = max(dhdl_avg)
min_y = min(dhdl_avg)

dhdl = np.zeros( (len(xvg_res), len(t)), dtype=float )
for k in range(len(xvg_sol)) :
    _, y = read_xvg_file( xvg_sol[k] )
    dhdl[k,:] = np.array(y)
dhdl_avg = np.mean(dhdl, axis=0)
delta_h = dt*np.cumsum(dhdl_avg)
diff_H_sol = np.mean(dhdl_avg)
print("diff_H_sol = "+str(diff_H_sol)+" kJ/mol")

diff_H_ws = diff_H_sol - diff_H_res - diff_H_sur

print("excess WS  = "+str(diff_H_ws)+" kJ/mol")

ax2.plot(lam, dhdl_avg, 'c-', linewidth=0.25)
ax2.plot(lam, delta_h, 'r-', linewidth=2.5)
ax2.plot(lam, np.zeros(len(lam)), 'k--', linewidth=2.5)

ax2.set_title('Solvated')
ax2.set_xlabel(r'$\lambda$')

max_y = max(max_y,max(dhdl_avg))
min_y = min(min_y,min(dhdl_avg))

ax1.set_ylim(min_y, max_y)
ax2.set_ylim(min_y, max_y)

plt.show()
