import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import mdconf as md

import os

# Times [ps]
Tinit = 1000.0
Tfin  = 2000.0
dt    = 12.5
n0    = int(Tinit/dt) 

folder = "/home/michele/python_for_md/TemperatureMap"
file_name = folder+"/conf_large.gro"

N = int((Tfin-Tinit)/dt)
N_sub = 2*8520
half_plane = 0.5*30.36600

oxygen_posx = np.zeros(N_sub)
oxygen_posy = np.zeros(N_sub)
oxygen_posz = np.zeros(N_sub)
oxygen_vel2 = np.zeros(N_sub)

for n in range(N) :
    
    # Dumping molecular snapshot
    t = Tinit + n*dt
    trjconv_command = "echo 2 | gmx trjconv -dump "+str(int(t)) \
            +" -f "+folder+"/traj.trr -s "+folder+"/system.tpr -o "+file_name
    os.system(trjconv_command)
    
    # Init list of position/velocities
    opx_tmp = []
    opy_tmp = []
    opz_tmp = []
    ov2_tmp = []
    
    # Reading snapshot
    n_lines = md.count_line( file_name )
    input_file = open(file_name, 'r')
    nn = 0
    with input_file as f :
        for line in f :
            nn += 1
            if nn > 2 and nn < n_lines :
                line_data = md.read_gro_line(line)
                if line_data[2] == "O1" or line_data[2] == "O2" :
                    opx_tmp.append(line_data[4])
                    opy_tmp.append(line_data[5])
                    opz_tmp.append(line_data[6])
                    ov2_tmp.append(line_data[7]**2+line_data[8]**2+line_data[9]**2)
    input_file.close()

    # Aggregating results
    oxygen_posx += np.array(opx_tmp)
    oxygen_posy += np.array(opy_tmp)
    oxygen_posz += np.array(opz_tmp)
    oxygen_vel2 += np.array(ov2_tmp)

    # Cleaning molecular sampshot
    os.system("rm "+folder+"/#*")

oxygen_posx /= N
oxygen_posy /= N
oxygen_posz /= N
oxygen_vel2 /= N

plt.scatter(oxygen_posx, oxygen_posz, c=oxygen_vel2, cmap=cm.jet)
plt.colorbar()
plt.title("Sub. oxygen $<v^2>$", fontsize=20.0)
plt.xticks(fontsize=15.0)
plt.yticks(fontsize=15.0)
plt.xlabel("x [nm]", fontsize=15.0)
plt.ylabel("z [nm]", fontsize=15.0)
plt.show()
