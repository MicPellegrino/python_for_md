import mdconf as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

file_name0 = '/home/michele/BeskowDiag/ShearSubstrate/104deg_NO_walltreat/confout_3ns.gro'
file_name1 = '/home/michele/BeskowDiag/ShearSubstrate/104deg_YES_walltreat/confout_3ns.gro'

N = []
N.append( md.count_line(file_name0) )
N.append( md.count_line(file_name1) )

# From inpput file
half_plane = 0.5*30.36600
d_so = 0.151
d_oo = 2.0*d_so

input_file = []
input_file.append( open(file_name0, 'r') )
input_file.append( open(file_name1, 'r') )

fig, axes = plt.subplots(nrows=2, ncols=1)

for idx in range(2) :
    
    N_lin = N[idx]

    n = 0
    silica_x = []
    silica_y = []
    silica_z = []
    oxigen1_x = []
    oxigen1_y = []
    oxigen1_z = []
    oxigen2_x = []
    oxigen2_y = []
    oxigen2_z = []
    
    for line in input_file[idx] :
        if n>=2 and n<N_lin-1 :
            line_data = md.read_gro_line(line)
            if line_data[6]<half_plane :
                if line_data[2] == "SI" :
                    silica_x.append(line_data[4])
                    silica_y.append(line_data[5])
                    silica_z.append(line_data[6])
                if line_data[2] == "O1" :
                    oxigen1_x.append(line_data[4])
                    oxigen1_y.append(line_data[5])
                    oxigen1_z.append(line_data[6])
                if line_data[2] == "O2" :
                    oxigen2_x.append(line_data[4])
                    oxigen2_y.append(line_data[5])
                    oxigen2_z.append(line_data[6])
        n += 1

    input_file[idx].close()

    silica_x = np.array(silica_x)
    silica_y = np.array(silica_y)
    silica_z = np.array(silica_z)
    oxigen1_x = np.array(oxigen1_x)
    oxigen1_y = np.array(oxigen1_y)
    oxigen1_z = np.array(oxigen1_z)
    oxigen2_x = np.array(oxigen2_x)
    oxigen2_y = np.array(oxigen2_y) 
    oxigen2_z = np.array(oxigen2_z)

    compression = d_oo - np.sqrt( (oxigen1_x-oxigen2_x)**2 + (oxigen1_y-oxigen2_y)**2 + (oxigen1_z-oxigen2_z)**2 )
    
    ax = axes.flat[idx]
    im = ax.scatter(silica_x, silica_y, c=compression, vmin=0.000, vmax=0.002, cmap=cm.jet, s=5.0)
    ax.set_xlim([min(silica_x), max(silica_x)])
    ax.set_ylim([min(silica_y), max(silica_y)])
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x [nm]', fontsize=15.0)
    ax.set_ylabel('y [nm]', fontsize=15.0)
    if idx == 0 :
        ax.set_title('Compress [nm], wall interactions', fontsize=15.0)
    else :
        ax.set_title('Compress [nm], no wall interactions', fontsize=15.0)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)

plt.show()

"""
plt.scatter(silica_x, silica_y, c=theta, cmap=cm.jet, s=5.0)
plt.colorbar()
plt.xlim([min(silica_x), max(silica_x)])
plt.ylim([min(silica_y), max(silica_y)])
plt.title('Deviation from vertical position [deg]', fontsize=20.0)
plt.xticks(fontsize=15.0)
plt.yticks(fontsize=15.0)
plt.xlabel('x [nm]', fontsize=15.0)
plt.ylabel('y [nm]', fontsize=15.0)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
"""
