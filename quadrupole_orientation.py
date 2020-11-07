import mdconf as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

file_name = 'conf_q2_walls.gro'

N_lin = md.count_line(file_name)

# From inpput file
half_plane = 0.5*30.36600;

input_file = open(file_name, 'r')

n = 0
silica_x = []
silica_y = []
oxigen1_x = []
oxigen1_y = []
oxigen1_z = []
oxigen2_x = []
oxigen2_y = []
oxigen2_z = []
for line in input_file :
    if n>=2 and n<N_lin-1 :
        line_data = md.read_gro_line(line)
        if line_data[6]<half_plane :
            if line_data[2] == "SI" :
                silica_x.append(line_data[4])
                silica_y.append(line_data[5])
            if line_data[2] == "O1" :
                oxigen1_x.append(line_data[4])
                oxigen1_y.append(line_data[5])
                oxigen1_z.append(line_data[6])
            if line_data[2] == "O2" :
                oxigen2_x.append(line_data[4])
                oxigen2_y.append(line_data[5])
                oxigen2_z.append(line_data[6])
    n += 1

input_file.close()

silica_x = np.array(silica_x)
silica_y = np.array(silica_y)
oxigen1_x = np.array(oxigen1_x)
oxigen1_y = np.array(oxigen1_y)
oxigen1_z = np.array(oxigen1_z)
oxigen2_x = np.array(oxigen2_x)
oxigen2_y = np.array(oxigen2_y) 
oxigen2_z = np.array(oxigen2_z)

diff_z = oxigen1_z-oxigen2_z
norm = np.sqrt( (oxigen1_x-oxigen2_x)**2 + (oxigen1_y-oxigen2_y)**2 + (oxigen1_z-oxigen2_z)**2 )
theta = np.arccos( diff_z/norm )

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
