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

label = ['00', '01', '02', '03', '04', '05', '06']
a = np.linspace(0.0, 0.7, 8)

xvg_res = []
xvg_sol = []
gro_init_res = []
gro_fin_res = []
for l in label :
    xvg_res.append('/home/michele/python_for_md/FreeEnergyCorrugate/FE'+l+'_res/dhdl.xvg')
    xvg_sol.append('/home/michele/python_for_md/FreeEnergyCorrugate/FE'+l+'_sol/dhdl.xvg')
    gro_init_res.append('/home/michele/python_for_md/FreeEnergyCorrugate/FE'+l+'_res/init_conf.gro')
    gro_fin_res.append('/home/michele/python_for_md/FreeEnergyCorrugate/FE'+l+'_res/confout.gro')

conf_init = md.Configuration()
conf_init.input(gro_init_res[0])
conf_fin = md.Configuration()
conf_fin.input(gro_fin_res[0])
Ly = (conf_init.box_yy)*1e-9
Lx0 = (conf_init.box_xx)*1e-9
Lx1 = (conf_fin.box_xx)*1e-9
diff_H_sur = [ (2.0*N_A*gamma*Ly*(Lx1-Lx0))*1e-3, ]

theta_0 = 37.8
delta_F_wet_dry = -2.0*N_A*gamma*np.cos(np.deg2rad(theta_0))*Ly*Lx0*1e-3

k_example = 4

for k in range(1, len(label)) :
    conf_init.input(gro_init_res[k])
    conf_fin.input(gro_fin_res[k])
    Lx0 = (conf_init.box_xx)*1e-9
    Lx1 = (conf_fin.box_xx)*1e-9
    diff_H_sur.append( (2.0*N_A*gamma*Ly*(Lx1-Lx0))*1e-3 )

t, _ = read_xvg_file( xvg_res[0] )
lam = np.array(t) / max(t)
dl = lam[1]-lam[0]

fig1, (ax1, ax2) = plt.subplots(1, 2)

diff_H_res = []
for k in range(len(label)) :
    _, dhdl_avg = read_xvg_file( xvg_res[k] )
    delta_h = dl*np.cumsum(dhdl_avg)
    diff_H_res.append(delta_h[-1])
    if k == k_example :
        ax1.plot(t, dhdl_avg, 'g-', linewidth=0.25)
        ax1.plot(t, np.zeros(len(t)), 'k--', linewidth=2.5)
        ax1.plot(t, delta_h, 'r-', linewidth=2.5, label=r'integral over $d\lambda$')
        ax1.set_title('Dry', fontsize=35.0)
        ax1.set_xlabel('$t$ [ps]', fontsize=25.0)
        ax1.set_ylabel(r'$<dH/d\lambda>$ [kJ/mol]', fontsize=25.0)
        ax1.legend(fontsize=25.0)
        ax1.tick_params(axis='x', labelsize=20.0)
        ax1.tick_params(axis='y', labelsize=20.0)
        y_min = min(dhdl_avg)
        y_max = max(dhdl_avg)

diff_H_sol = []
for k in range(len(label)) :
    _, dhdl_avg = read_xvg_file( xvg_sol[k] )
    delta_h = dl*np.cumsum(dhdl_avg)
    diff_H_sol.append(delta_h[-1])
    if k == k_example :
        ax2.plot(t, dhdl_avg, 'c-', linewidth=0.25)
        ax2.plot(t, np.zeros(len(t)), 'k--', linewidth=2.5)
        ax2.plot(t, delta_h, 'r-', linewidth=2.5)
        ax2.set_title('Wet', fontsize=35.0)
        ax2.set_xlabel('$t$ [ps]', fontsize=25.0)
        ax2.tick_params(axis='x', labelsize=20.0)
        ax2.tick_params(axis='y', labelsize=20.0)
        y_min = min( y_min, min(dhdl_avg))
        y_max = max( y_max, max(dhdl_avg))
        ax1.set_xlim([0.0, max(t)])
        ax1.set_ylim([y_min, y_max])
        ax2.set_xlim([0.0, max(t)])
        ax2.set_ylim([y_min, y_max])

plt.show()

diff_H_sol.insert(0, 0.0)
diff_H_res.insert(0, 0.0)
diff_H_sur.insert(0, 0.0)

diff_H_sol = np.array(diff_H_sol)
diff_H_res = np.array(diff_H_res)
diff_H_sur = np.array(diff_H_sur)

diff_H_ws = diff_H_sol - diff_H_res - diff_H_sur

plt.plot(a, diff_H_sol, 'ko-.', linewidth=1.5, markersize=10.0, label='total')
plt.plot(a, diff_H_res, 'gs-.', linewidth=1.5, markersize=10.0, label='pos. res.')
plt.plot(a, diff_H_sur, 'bx-.', linewidth=1.5, markersize=10.0, label='liquid/vapour')
plt.plot(a, diff_H_ws,  'rD-.', linewidth=3.0, markersize=15.0, label='excess')
plt.plot([0.0, 0.75], [delta_F_wet_dry, delta_F_wet_dry], 'm--', linewidth=3.0, label=r'$\delta$ wet-dry')
plt.title('Free energy analysis', fontsize=35.0)
plt.legend(fontsize=25.0)
plt.xlim(0.0, 0.75)
plt.xticks(fontsize=20.0)
plt.yticks(fontsize=20.0)
plt.ylabel('$\Delta F$ [kJ/mol]', fontsize=25.0)
plt.xlabel('$a=hk$ [-1]', fontsize=25.0)
plt.show()
