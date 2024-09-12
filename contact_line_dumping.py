import numpy as np
import os

input_folder = '/home/michele/BeskowDiag/Select_Q4_C0015'
output_folder = '/home/michele/python_for_md/ContactLineTracingQ4/RECL'

ndt = 5

# traj_file = input_folder+'/traj_22ns_26ns.trr'
traj_file = input_folder+'/trj_recl_52_56.gro'
system_file = input_folder+'/system.tpr'
index_file = input_folder+'/recl_54.ndx'

"""
T0 = 22000+1130
deltaT = 4000-1130
N = int(deltaT/ndt)
for i in range(N) :
    conf_file = output_folder+"/"+str(i*ndt).zfill(4)+"ps.gro"
    cl = "gmx trjconv -f "+traj_file+" -s "+system_file+" -n "+index_file+" -b "+str(T0+(i-1)*ndt)+" -dump "+str(T0+i*ndt)+" -o "+conf_file
    os.system(cl)
"""

n_atoms = 894
T = 4000
dt = 5
N = int(T/dt)
fin = open(traj_file, 'r')
n = 0
i = 0
for line in fin :
    if n == (n_atoms+3) :
        n = 0
        i += 1
        fout.close()
    if n == 0 :
        fout = open(output_folder+'/'+str(int(dt*i)).zfill(4)+'ps.gro', 'w')
    fout.write(line)
    n += 1
fin.close()
