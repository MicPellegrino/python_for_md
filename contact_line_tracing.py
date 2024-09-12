import numpy as np
import mdconf_oop as md
import matplotlib.pyplot as plt

t_max = 4000
dt = 5
m = int(t_max/dt)
cmap=plt.get_cmap("coolwarm")
aa = 0.35
ms = 7.5
lw = 1.75
fs = 30.0
folder_root = '/home/michele/python_for_md/ContactLineTracingQ2'

input_folder = folder_root+'/RECR'
cl = md.Configuration()
cl.input(input_folder+'/0000ps.gro')
cl.print_info()
n = 0
for a in cl.res_atomlist['SOL'] :
    if a.label=='OW' : 
        n += 1
print("We found "+str(n)+" oxygen atoms, let's follow them")
oxygen_px = np.zeros(n, dtype=float)
oxygen_pz = np.zeros(n, dtype=float)
com_x = []
com_z = []
for i in range(m) :
    cl = md.Configuration()
    cl.input(input_folder+'/'+str(i*dt).zfill(4)+'ps.gro')
    k = 0
    for a in cl.res_atomlist['SOL'] :
        if a.label=='OW' : 
            oxygen_px[k] = a.pos_x
            oxygen_pz[k] = a.pos_z
            k += 1
    com_x.append(np.mean(oxygen_px))
    com_z.append(np.mean(oxygen_pz))
    plt.plot(oxygen_px, oxygen_pz, '.', color=cmap(i/m), alpha=aa, markeredgecolor='None', markersize=ms)
    oxygen_px = np.zeros(n, dtype=float)
    oxygen_pz = np.zeros(n, dtype=float)
plt.plot(com_x, com_z, 'k-', linewidth=lw)

input_folder = folder_root+'/RECL'
cl = md.Configuration()
cl.input(input_folder+'/0000ps.gro')
cl.print_info()
n = 0
for a in cl.res_atomlist['SOL'] :
    if a.label=='OW' : 
        n += 1
print("We found "+str(n)+" oxygen atoms, let's follow them")
oxygen_px = np.zeros(n, dtype=float)
oxygen_pz = np.zeros(n, dtype=float)
com_x = []
com_z = []
for i in range(m) :
    cl = md.Configuration()
    cl.input(input_folder+'/'+str(i*dt).zfill(4)+'ps.gro')
    k = 0
    for a in cl.res_atomlist['SOL'] :
        if a.label=='OW' : 
            oxygen_px[k] = a.pos_x
            oxygen_pz[k] = a.pos_z
            k += 1
    com_x.append(np.mean(oxygen_px))
    com_z.append(np.mean(oxygen_pz))
    plt.plot(oxygen_px, oxygen_pz, '.', color=cmap(i/m), alpha=aa, markeredgecolor='None', markersize=ms)
    oxygen_px = np.zeros(n, dtype=float)
    oxygen_pz = np.zeros(n, dtype=float)
plt.plot(com_x, com_z, 'k-', linewidth=lw)

input_folder = folder_root+'/ADVR'
cl = md.Configuration()
cl.input(input_folder+'/0000ps.gro')
cl.print_info()
n = 0
for a in cl.res_atomlist['SOL'] :
    if a.label=='OW' : 
        n += 1
print("We found "+str(n)+" oxygen atoms, let's follow them")
oxygen_px = np.zeros(n, dtype=float)
oxygen_pz = np.zeros(n, dtype=float)
com_x = []
com_z = []
for i in range(m) :
    cl = md.Configuration()
    cl.input(input_folder+'/'+str(i*dt).zfill(4)+'ps.gro')
    k = 0
    for a in cl.res_atomlist['SOL'] :
        if a.label=='OW' : 
            oxygen_px[k] = a.pos_x
            oxygen_pz[k] = a.pos_z
            k += 1
    com_x.append(np.mean(oxygen_px))
    com_z.append(np.mean(oxygen_pz))
    plt.plot(oxygen_px, oxygen_pz, '.', color=cmap(i/m), alpha=aa, markeredgecolor='None', markersize=ms)
    oxygen_px = np.zeros(n, dtype=float)
    oxygen_pz = np.zeros(n, dtype=float)
plt.plot(com_x, com_z, 'k-', linewidth=lw)

input_folder = folder_root+'/ADVL'
cl = md.Configuration()
cl.input(input_folder+'/0000ps.gro')
cl.print_info()
n = 0
for a in cl.res_atomlist['SOL'] :
    if a.label=='OW' : 
        n += 1
print("We found "+str(n)+" oxygen atoms, let's follow them")
oxygen_px = np.zeros(n, dtype=float)
oxygen_pz = np.zeros(n, dtype=float)
com_x = []
com_z = []
for i in range(m) :
    cl = md.Configuration()
    cl.input(input_folder+'/'+str(i*dt).zfill(4)+'ps.gro')
    k = 0
    for a in cl.res_atomlist['SOL'] :
        if a.label=='OW' : 
            oxygen_px[k] = a.pos_x
            oxygen_pz[k] = a.pos_z
            k += 1
    com_x.append(np.mean(oxygen_px))
    com_z.append(np.mean(oxygen_pz))
    plt.plot(oxygen_px, oxygen_pz, '.', color=cmap(i/m), alpha=aa, markeredgecolor='None', markersize=ms)
    oxygen_px = np.zeros(n, dtype=float)
    oxygen_pz = np.zeros(n, dtype=float)
plt.plot(com_x, com_z, 'k-', linewidth=lw, label='local c.o.m.')

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-2.0, vmax=2.0))
time_scale = plt.colorbar(sm)
time_scale.set_label('time [ns]', rotation=270, fontsize=0.75*fs)
time_scale.ax.tick_params(labelsize=0.5*fs)
plt.title('Molecular transport to contact lines', fontsize=fs)
plt.xlabel('x [nm]', fontsize=0.75*fs)
plt.ylabel('y [nm]', fontsize=0.75*fs)
plt.axis('square')
plt.legend(fontsize=0.75*fs)
plt.xticks(fontsize=0.5*fs)
plt.yticks(fontsize=0.5*fs)
plt.xlim([30.0, 130.0])
plt.ylim([0.0,  30.63400])
plt.show()
