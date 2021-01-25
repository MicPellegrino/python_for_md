import numpy as np
import os

import mdconf as md

# Nanometers
Lx = 50.0000
Ly = 50.0000
Lz = 30.6340

N_h2o_target = 174

dist = 10.0*Lz - 2.0*4.50*np.sqrt(2.0/3.0)
ni = int( (10.0*Lx)/4.50 )
nj = int( (10.0*Ly)/(4.50*0.5*np.sqrt(3.0)) )
nk = 1

folder = '/home/michele/python_for_md'
output = folder+'/vapour_layers.pdb'

md.quad_double_layer(dist, ni, nj, nk, output, sep = 4.990000)

command = "gmx editconf -c -f "+output+" -o "+folder+"/vapour_layers.gro"
print(command)
os.system(command)

solvate = folder+"/vapour_solvate.gro"
command = "gmx solvate -cp "+folder+"/vapour_layers.gro -cs "+folder+"/wat_equil.gro -radius 4.0 -scale 8.0 -o "+solvate
print(command)
os.system(command)

##########################################################
# Script imposing the correct number of molecules inside #
##########################################################

init_conf = folder+"/vapour_init.gro"

n_lines = md.count_line( solvate )

f_in = open( solvate, 'r' )
idx = 0
n_h2o_atoms = 0 
for line in f_in :
    idx += 1
    if idx > 2 and idx < n_lines :
        line_data = md.read_gro_line( line )
        if line_data[1] == "SOL" :
            n_h2o_atoms += 1
f_in.close()

N_h2o_0 = int(n_h2o_atoms/3)
r = N_h2o_0/N_h2o_target

n_atoms_mol = 3
mol_line = ['', '', '']
f_in = open( solvate, 'r' )
f_out = open( init_conf, 'w' )
idx = 0
n_atoms = 0
mol_count = 0
header = ''
n = 0
for line in f_in :
    idx += 1
    if idx > 2 and idx < n_lines :
        line_data = md.read_gro_line( line )
        if line_data[1] == "SOL" :
            if (n//n_atoms_mol)%int(r) == 0 :
                f_out.write(line)
                n_atoms += 1
            n += 1
        else :
            f_out.write(line)
            n_atoms += 1
    elif idx == n_lines :
        f_out.write("%10.5f%10.5f%10.5f\n" % ( Lx, Ly, Lz ) )
    elif idx == 1 :
        header = line
        f_out.write(line)
    else :
        f_out.write(line)
f_out.close()
f_in.close()
f_out = open( init_conf, 'r+')
f_out.write(header)
f_out.write(" "+str(n_atoms)+"\t")
f_out.close()
