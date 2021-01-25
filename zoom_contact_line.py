import os
import mdconf as md

# Times [ps]
Tinit = 6250.0
Tfin  = 7500.0
dt    = 12.5
n0    = int(Tinit/dt) 

# Window [nm]
x0 = 87.06375
y0 = 0.00
z0 = 0.00
x1 = 100.0434375
y1 = 2.5
z1 = 5.993289473684211

folder = "/home/michele/BeskowDiag/Shear123deg"
N = int((Tfin-Tinit)/dt)
for n in range(N) :
    t = Tinit + n*dt
    trjconv_command = "echo 0 | gmx trjconv -dump "+str(int(t))+" -f "+folder+"/traj.trr -s "+folder+"/system.tpr -o "+folder+"/conf_large.gro"
    os.system(trjconv_command)
    md.crop_xyz( x0, y0, z0, x1, y1, z1, folder+"/conf_large.gro", folder+"/zoom.gro")
    editconf_command = "gmx editconf -f "+folder+"/zoom.gro -o "+folder+"/ConfSnapshots/zoom_"+str(n0+n).zfill(5)+".pdb"
    os.system(editconf_command)
    os.system("rm "+folder+"/#*")
