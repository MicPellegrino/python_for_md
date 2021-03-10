import os
import mdconf as md

# Times [ps]
Tinit = 4000.0
Tfin  = 5000.0
dt    = 20.0
n0    = int(Tinit/dt) 

# Window [nm]
x0 = 50.00
y0 = 0.00
z0 = 10.00
x1 = 110.0
y1 = 4.67650
z1 = 10.50

folder = "/home/michele/BeskowDiag/VelocityDistribution"
N = int((Tfin-Tinit)/dt)
for n in range(N) :
    t = Tinit + n*dt
    trjconv_command = "echo 4 | gmx trjconv -dump "+str(int(t))+" -f "+folder+"/traj.trr -s "+folder+"/system.tpr -o "+folder+"/conf_large.gro"
    os.system(trjconv_command)
    md.crop_xyz( x0, y0, z0, x1, y1, z1, folder+"/conf_large.gro", folder+"/ConfSnapshots/zoom_"+str(n0+n).zfill(5)+".gro")
    os.system("rm "+folder+"/#*")
