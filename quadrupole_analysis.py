import mdconf as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

label_vec = ['2', '3', '4', '5']
# label_vec = ['5']

folder_name0 = '/home/michele/BeskowDiag/ShearSubstrate/104deg_NO_walltreat/'
folder_name1 = '/home/michele/BeskowDiag/ShearSubstrate/104deg_YES_walltreat/'

N_sub = 4260
half_plane = 0.5*30.36600
d_oo = 0.151
z_restrain = 0.366

fig, axes = plt.subplots(nrows=2, ncols=1)

for idx in range(2) :
    
    silica_x  = np.zeros(N_sub)
    silica_y  = np.zeros(N_sub)
    silica_z  = np.zeros(N_sub)
    oxigen1_x = np.zeros(N_sub)
    oxigen1_y = np.zeros(N_sub)
    oxigen1_z = np.zeros(N_sub)
    oxigen2_x = np.zeros(N_sub)
    oxigen2_y = np.zeros(N_sub)
    oxigen2_z = np.zeros(N_sub)

    for l in label_vec :
        
        file_name = ''
        if idx == 0 :
            file_name = folder_name0+'confout_'+l+'ns.gro'
        else :
            file_name = folder_name1+'confout_'+l+'ns.gro'

        N_lin = md.count_line(file_name)
        
        input_file = open(file_name, 'r')

        n = 0
        s_x = []
        s_y = []
        s_z = []
        o1_x = []
        o1_y = []
        o1_z = []
        o2_x = []
        o2_y = []
        o2_z = []
    
        for line in input_file :
            if n>=2 and n<N_lin-1 :
                line_data = md.read_gro_line(line)
                if line_data[6]<half_plane :
                    if line_data[2] == "SI" :
                        s_x.append(line_data[4])
                        s_y.append(line_data[5])
                        s_z.append(line_data[6])
                    if line_data[2] == "O1" :
                        o1_x.append(line_data[4])
                        o1_y.append(line_data[5])
                        o1_z.append(line_data[6])
                    if line_data[2] == "O2" :
                        o2_x.append(line_data[4])
                        o2_y.append(line_data[5])
                        o2_z.append(line_data[6])
            n += 1

        input_file.close()

        s_x = np.array(s_x)
        s_y = np.array(s_y)
        s_z = np.array(s_z)
        o1_x = np.array(o1_x)
        o1_y = np.array(o1_y)
        o1_z = np.array(o1_z)
        o2_x = np.array(o2_x)
        o2_y = np.array(o2_y) 
        o2_z = np.array(o2_z)

    silica_x  += s_x / len(label_vec)
    silica_y  += s_y / len(label_vec)
    silica_z  += s_z / len(label_vec)
    oxigen1_x += o1_x / len(label_vec)
    oxigen1_y += o1_y / len(label_vec)
    oxigen1_z += o1_z / len(label_vec)
    oxigen2_x += o2_x / len(label_vec)
    oxigen2_y += o2_y / len(label_vec)
    oxigen2_z += o2_z / len(label_vec)

    diff_z = oxigen1_z-oxigen2_z
    norm = np.sqrt( (oxigen1_x-oxigen2_x)**2 + (oxigen1_y-oxigen2_y)**2 + (oxigen1_z-oxigen2_z)**2 )
    theta = np.arccos( diff_z/norm )

    compression = d_oo - np.sqrt( (oxigen1_x-oxigen2_x)**2 + (oxigen1_y-oxigen2_y)**2 + (oxigen1_z-oxigen2_z)**2 )
    shift_z = silica_z - z_restrain

    ax = axes.flat[idx]
    
    im = ax.scatter(silica_x, silica_y, c=shift_z, vmin=0.000, vmax=0.175, cmap=cm.jet, s=5.0)
    if idx == 0 :
        ax.set_title('Shift [nm], wall interactions', fontsize=15.0)
    else :
        ax.set_title('Shift [nm], no wall interactions', fontsize=15.0)

    """
    im = ax.scatter(silica_x, silica_y, c=compression, vmin=0.000, vmax=0.10, cmap=cm.jet, s=5.0)
    if idx == 0 :
        ax.set_title('Compress [nm], wall interactions', fontsize=15.0)
    else :
        ax.set_title('Compress [nm], no wall interactions', fontsize=15.0)
    """

    """
    im = ax.scatter(silica_x, silica_y, c=theta, vmin=0.0, vmax=1.5 ,cmap=cm.jet, s=5.0)
    if idx == 0 :
        ax.set_title('Orient [deg], wall interactions', fontsize=15.0)
    else :
        ax.set_title('Orient [deg], no wall interactions', fontsize=15.0)
    """

    ax.set_xlim([min(silica_x), max(silica_x)])
    ax.set_ylim([min(silica_y), max(silica_y)])
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x [nm]', fontsize=15.0)
    ax.set_ylabel('y [nm]', fontsize=15.0)

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
