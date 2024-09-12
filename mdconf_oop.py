import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import scipy as sc
import scipy.special

# DEFINITION THAT SHOULD NOT BE CHANGED BY THE USER #

# Silicon-oxygen distance in silica quadrupoles
d_so = 0.151
# Bond parameters for alkane grafting (carbon-carbon and carbon-hydrogen)
d_cc = 0.1529
d_ch = 0.1090
# Atom types in silica quadrupoles
asl = "SI"
ao1 = "O1"
ao2 = "O2"
# Atom name for carbon in grafted silica
ac1 = "CS"
# Atom types for LJ walls
ljw = "CUB"
# Hexagonal silica lattice parameters
alpha_1 = np.sqrt(3.0/4.0)
alpha_2 = np.sqrt(2.0/3.0)
# NUmber of steps to differnciate local geometry on rough substrates
inner_steps = 5
# Effective surface as function of roughness parameter
r_rough = lambda a2 : (2.0/np.pi) * np.sqrt(a2+1.0) * sc.special.ellipe(a2/(a2+1.0))

# Bonds and angles for amorphous silica (from CHARMM)
d_ct2_oh1 = 0.142
d_oh1_h = 0.096
a_ct2_oh1_h = 106.0000

#####################################################

"""
    Utilities
"""

def array_from_file( filename ):
    my_list = []
    with open(filename, 'r') as f:
        for line in f:
            my_list.append(float(line.split()[0]))
    return np.array(my_list)

def count_line( file_name ) :
    n_lines = 0
    f_in = open( file_name, 'r')
    for line in f_in :
        n_lines += 1
    f_in.close()
    return n_lines

def non_null_float( s ) :
    if s == '' or s == '\n' :
        return 0.0
    else :
        return float( s )

def read_gro_line( line, read_vel=True ) :
    line_data = [None]*(7+3*read_vel)
    line_data[0] = int(line[0:5].strip())                   # Residue serial
    line_data[1] = str(line[5:10].strip())                  # Residue name
    line_data[2] = str(line[10:15].strip())                 # Atom name
    line_data[3] = int(line[15:20].strip())                 # Atom serial
    line_data[4] = float(line[20:28].strip())               # Pos X
    line_data[5] = float(line[28:36].strip())               # Pos Y
    line_data[6] = float(line[36:44].strip())               # Pos Z
    if read_vel :
        line_data[7] = non_null_float( line[44:52].strip() )    # Vel X
        line_data[8] = non_null_float( line[52:60].strip() )    # Vel Y
        line_data[9] = non_null_float( line[60:68].strip() )    # Vel Z
    return line_data

def read_pdb_line( line ) :
    # Incipit = ATOM
    line_data = [None]*7
    line_data[0] = int(line[8:11].strip())      # Atom serial
    line_data[1] = str(line[12:16].strip())     # Atom name
    line_data[2] = str(line[17:20].strip())     # Residue name
    line_data[3] = int(line[22:26].strip())     # Residue serial
    line_data[4] = float(line[30:38].strip())   # Pos X
    line_data[5] = float(line[38:46].strip())   # Pos Y
    line_data[6] = float(line[46:54].strip())   # Pos Z
    return line_data

def read_lmp_pos( cols ) :
    line_data = [None]*7
    line_data[0] = int(cols[0])         # Atom number
    line_data[1] = str(cols[1])         # Atom group (unused)
    line_data[2] = str(cols[2])         # Atom type
    line_data[3] = float(cols[3])       # q (???)
    line_data[4] = 0.1*float(cols[4])       # Pos X
    line_data[5] = 0.1*float(cols[5])       # Pos Y
    line_data[6] = 0.1*float(cols[6])       # Pos Z
    return line_data

def read_lmp_vel( cols ) :
    line_data = [None]*4
    line_data[0] = int(cols[0])         # Atom number
    line_data[1] = 0.1*float(cols[1])       # Vel x
    line_data[2] = 0.1*float(cols[2])       # Vel y
    line_data[3] = 0.1*float(cols[3])       # Vel z
    return line_data


"""
    Classes
"""

class Atom :

    def __init__(self, l='', n=0, px=0.0, py=0.0, pz=0.0, vx=0.0, vy=0.0, vz=0.0):
        self.label = l
        self.index = n
        self.pos_x = px
        self.pos_y = py
        self.pos_z = pz
        self.vel_x = vx
        self.vel_y = vy
        self.vel_z = vz

    def gro_line(self, res_idx, res_name) :
        lg = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % ( res_idx, res_name, self.label, self.index % 100000,
                self.pos_x, self.pos_y, self.pos_z, self.vel_x, self.vel_y, self.vel_z )
        return lg

    def pdb_line(self, res_idx, res_name) :
        lp = "%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" % ("ATOM", self.index % 100000, self.label, res_name, res_idx, 
                self.pos_x, self.pos_y, self.pos_z)
        return lp


class Configuration :
    
    """
        Initializes an empty box
    """
    def __init__(self, Lx=0.0, Ly=0.0, Lz=0.0, h="") :
        
        print("Initializing an empty simulation box")
        
        # Instance-level memebers
        self.box_xx = Lx
        self.box_yy = Ly
        self.box_zz = Lz
        self.header = h
        self.n_atoms = 0
        self.header = ""

        # Dictionary containing list of atoms per each residue
        self.res_atomlist = dict()
        # Dictionary containing the set of atoms labels for each residue
        self.res_nomicon = dict()

    def pcb_dist(self, a1, a2) :

        dx = a2.pos_x - a1.pos_x
        dy = a2.pos_y - a1.pos_y
        dz = a2.pos_z - a1.pos_z

        if dx > 0.5*self.box_xx :
            dx -= self.box_xx
        elif dx <= -0.5*self.box_xx :
            dx += self.box_xx
        if dy > 0.5*self.box_yy :
            dy -= self.box_yy
        elif dy <= -0.5*self.box_yy :
            dy += self.box_yy
        if dz > 0.5*self.box_zz :
            dz -= self.box_zz
        elif dz <= -0.5*self.box_zz :
            dz -= self.box_zz

        return np.sqrt(dx*dx+dy*dy+dz*dz)


    """
        Re-initializes configuration via input file
        (Python OOP is stupid and you can't have nice things like constructor overloading ...)
    """
    def input(self, file_name) :

        ext = file_name.split('.')[-1]

        if ext == 'gro':
            print("Reading .gro configuration file")

            n_lines = count_line(file_name)
            n = 0 
            in_file = open(file_name, 'r')
            for line in in_file :
                n += 1
                if n == 1 :
                    self.header = line.strip('\n') 
                elif n == 2 :
                    self.n_atoms = int(line)
                elif n == n_lines :
                    cols = line.split()
                    self.box_xx = float(cols[0])
                    self.box_yy = float(cols[1])
                    self.box_zz = float(cols[2])
                else :
                    line_data = read_gro_line(line)
                    if not( line_data[1] in self.res_atomlist.keys() ) :
                        self.res_atomlist[line_data[1]] = []
                        self.res_nomicon[line_data[1]] = set()
                    self.res_atomlist[line_data[1]].append( Atom(line_data[2], line_data[3], \
                        line_data[4], line_data[5], line_data[6], \
                        line_data[7], line_data[8], line_data[9]) )
                    self.res_nomicon[line_data[1]].add(line_data[2])
            in_file.close()

        elif ext == 'pdb':
            print("Reading .pdb configuration file")
            in_file = open(file_name, 'r')
            atom_idx = 0
            for line in in_file :
                if line[0:6] == "ATOM  " :
                    atom_idx += 1
                    line_data = read_pdb_line(line)
                    if not( line_data[2] in self.res_atomlist.keys() ) :
                        self.res_atomlist[line_data[2]] = []
                        self.res_nomicon[line_data[2]] = set()
                    self.res_atomlist[line_data[2]].append( Atom(line_data[1], atom_idx, \
                          0.1*line_data[4], 0.1*line_data[5], 0.1*line_data[6], 0.0, 0.0, 0.0) )
                    self.res_nomicon[line_data[2]].add(line_data[1])
                elif line[0:6] == "CRYST1" :
                    # UNSAFE!
                    cols = line.split()
                    self.box_xx = 0.1*float(cols[1])
                    self.box_yy = 0.1*float(cols[2])
                    self.box_zz = 0.1*float(cols[3])
                elif line[0:6] == "HEADER" :
                    self.header = line.strip('\n') 
            in_file.close()
            for resname in self.res_atomlist.keys() :
                    self.n_atoms += len(self.res_atomlist[resname])

        else :
            print("Unsupported extension, please supply either .gro or .pdb")


    """
    def read_lmp_pos( cols ) :
    line_data = [None]*7
    line_data[0] = int(cols[0])         # Atom number
    line_data[1] = str(cols[1])         # Atom group (unused)
    line_data[2] = str(cols[2])         # Atom type
    line_data[3] = float(cols[3])       # q (???)
    line_data[4] = float(cols[4])       # Pos X
    line_data[5] = float(cols[5])       # Pos Y
    line_data[6] = float(cols[6])       # Pos Z
    return line_data
    def read_lmp_vel( cols ) :
    line_data = [None]*4
    line_data[0] = int(cols[0])         # Atom number
    line_data[1] = float(cols[1])       # Vel x
    line_data[2] = float(cols[2])       # Vel y
    line_data[3] = float(cols[3])       # Vel z
    return line_data
    """

    """
        Input from LAMMPS data files (because for fucking reasons there is no toll that does LAMMPS to .gro/.pdb automatically!)
    """
    def input_lammps(self, file_name, resname='SIL') :

        print("Reading LAMMPS configuration file")

        n = 0
        n_lines = count_line(file_name)
        na = n_lines
        nv = n_lines
        nvv = 0
        in_file = open(file_name, 'r')

        self.res_atomlist[resname] = []
        self.res_nomicon[resname] = set()

        for line in in_file :
            if line.strip('\n') :
                n += 1
                cols = line.split()
                if cols[-1] == 'atoms' :
                    self.n_atoms = int(cols[0])
                elif cols[-1] == 'xhi' :
                    self.box_xx = 0.1*float(cols[1])
                elif cols[-1] == 'yhi' :
                    self.box_yy = 0.1*float(cols[1])
                elif cols[-1] == 'zhi' :
                    self.box_zz = 0.1*float(cols[1])
                elif cols[0] == 'Atoms' :
                    na = n
                elif n>na and n<nv and cols[0] != 'Velocities' :
                    line_data = read_lmp_pos(cols)
                    self.res_atomlist[resname].append( Atom(line_data[2], line_data[0], \
                        line_data[4], line_data[5], line_data[6], \
                        0.0, 0.0, 0.0) )
                    self.res_nomicon[resname].add(line_data[2])
                elif cols[0] == 'Velocities' :
                    nv = n
                elif n>nv :
                    line_data = read_lmp_vel(cols)
                    self.res_atomlist[resname][nvv].vx = line_data[1]
                    self.res_atomlist[resname][nvv].vy = line_data[2]
                    self.res_atomlist[resname][nvv].vz = line_data[3]
                    nvv += 1

        in_file.close()

        # DOUBLE-CHECK
        assert (self.n_atoms==len(self.res_atomlist[resname])), "Inconsistent number of atoms!"

    
    """
        Input sphere (hoping it's not killed!)
    """
    def input_sphere(self, file_name, radius, cx, cy, cz, resname='SOL') :

        print("Reading no. lines")
        n_lines = count_line(file_name)
        n = 0
        cut2 = radius**2

        print("Opening file")
        in_file = open(file_name, 'r')
        n_count = 0

        print("Initializing the structure")
        self.res_atomlist[resname] = []
        self.res_nomicon[resname] = set()
        loc_atom_list = [None] * 3

        for line in in_file :
           
            n += 1
            
            if n == 1 :
                self.header = line.strip('\n')
            
            elif n == n_lines :
                cols = line.split()
                self.box_xx = float(cols[0])
                self.box_yy = float(cols[1])
                self.box_zz = float(cols[2])
            
            elif n > 2 :
                
                line_data = read_gro_line(line)
                
                a = Atom(line_data[2], line_data[3], \
                    line_data[4], line_data[5], line_data[6], \
                    line_data[7], line_data[8], line_data[9])
                
                x = a.pos_x
                y = a.pos_y
                z = a.pos_z

                loc_atom_list[n_count] = a
                
                if n_count == 0 :
                    r2 = (x-cx)**2 + (y-cy)**2 + (z-cz)**2
                else :
                    r2 = max(r2, (x-cx)**2 + (y-cy)**2 + (z-cz)**2)

                n_count += 1
                if n_count == 3 :
                    n_count = 0

                    if r2 < cut2 :
                        
                        for k in range(3) :
                            
                            assert (loc_atom_list[k]!=None), "Colder than absolute zero!"
                            
                            self.res_atomlist[resname].append( loc_atom_list[k] )
                            self.res_nomicon[resname].add(line_data[2])
                            
                            if len(self.res_atomlist[resname])%10000 == 0 :
                                print("We just added another 10k atoms, yippe!")

        in_file.close()

        self.n_atoms += len(self.res_atomlist[resname])


    
    """
        Outputs as .gro configuration file
    """
    def output(self, file_name) :
        
        ext = file_name.split('.')[-1]

        if ext == 'gro':
            print("Writing .gro configuration file")
            out_file = open(file_name, 'w')
            out_file.write(self.header+'\n')
            out_file.write(str(self.n_atoms)+'\n')
            res_idx = 0
            for res_name in self.res_atomlist.keys() :
                res_idx += 1
                for atom in  self.res_atomlist[res_name] :
                    out_file.write(atom.gro_line(res_idx, res_name))
            out_file.write("%10.5f%10.5f%10.5f\n" % \
                ( self.box_xx, self.box_yy, self.box_zz ) )
            out_file.close()

        else :
            print("Only .gro configuration files are supported as output")

    """
        Simple output function for amorphous silica annealing simulations
    """
    def output_lammps(self, file_name, mass=None) :

        out_file = open(file_name, 'w')

        out_file.write(self.header+'\n')
        out_file.write('\n')

        out_file.write(str(self.n_atoms)+' atoms\n')
        nat = 0
        for resname in self.res_nomicon.keys() :
            nat += len(self.res_nomicon[resname])
        out_file.write(str(nat)+' atom types\n')
        out_file.write('\n')

        out_file.write("0.000000 "+'{:.6f}'.format(10*self.box_xx)+' xlo xhi\n')
        out_file.write("0.000000 "+'{:.6f}'.format(10*self.box_yy)+' ylo yhi\n')
        out_file.write("0.000000 "+'{:.6f}'.format(10*self.box_zz)+' zlo zhi\n')
        out_file.write('\n')

        out_file.write(' Masses\n')
        out_file.write('\n')
        if mass == None :
            for resname in self.res_nomicon.keys() :
                atl = sorted(self.res_nomicon[resname])
                for at in atl :
                    out_file.write(at+" 0.000000\n")
        else :
            for resname in self.res_nomicon.keys() :
                atl = sorted(self.res_nomicon[resname])
                for at in atl :
                    out_file.write(at+' {:.4f}'.format(mass[at])+" \n")
        out_file.write('\n')

        out_file.write('Atoms # full\n')
        out_file.write('\n')
        ri = 0
        na = 0
        for rn in self.res_atomlist.keys() :
            ri += 1
            for a in  self.res_atomlist[rn] :
                na += 1
                px = '{:.6f}'.format(10*a.pos_x)
                py = '{:.6f}'.format(10*a.pos_y)
                pz = '{:.6f}'.format(10*a.pos_z)
                out_file.write(str(na)+' '+str(ri)+' '+str(a.label)+' '+'0.000000 '+px+' '+py+' '+pz+' # '+str(a.label)+' '+rn+' \n')

        out_file.close()


    """
        Creates a flat LJ wall of nk layers
    """
    def lj_wall(self, z_wall, nk=3, sp=0.27, resname='SUB', cut_box=False) :

        dx = sp
        dy = sp*alpha_1
        dz = sp*alpha_2

        ni = int(self.box_xx/dx)
        nj = int(self.box_yy/dy)

        if cut_box :
            self.box_xx = ni*dx
            # Make sure that when cutting long y the number of layes is even (to avoid period images inconst.)
            nj += (nj%2)
            self.box_yy = nj*dy

        if not( resname in self.res_atomlist.keys() ) :
            self.res_atomlist[resname] = []
            self.res_nomicon[resname] = set()
            self.res_nomicon[resname].add(ljw)

        dx_y = dx/2
        dx_z = dx/2
        dy_z = dy/3
        x0 = dx/4
        y0 = dy/6
        z0 = z_wall-0.5*dy*nk

        wrap = lambda t, bt : t-bt*np.floor(t/bt)

        n = len(self.res_atomlist[resname])
        for k in range(nk):
            a = "CUB"
            z = z0 + k*dz
            for i in range(ni):
                for j in range(nj):
                    n += 1
                    y = y0 + j*dy + k*dy_z
                    x = x0 + i*dx + j*dx_y + k*dx_z
                    self.res_atomlist[resname].append( Atom(ljw, n, wrap(x,self.box_xx), wrap(y,self.box_yy), z, 0.0, 0.0, 0.0) )
                    self.n_atoms += 1

    """
        Creates a flat silica monolayer at the prescribed z coordinate
    """
    def silica_monolayer(self, z_wall, sp=0.450, resname='SUB', cut_box=False) :
        
        dx = sp
        dy = sp*alpha_1
        dx_y = dx/2.0
        
        ni = int(self.box_xx/dx)
        nj = int(self.box_yy/dy)

        if cut_box :
            self.box_xx = ni*dx
            # Make sure that when cutting long y the number of layes is even (to avoid period images inconst.)
            nj += (nj%2)
            self.box_yy = nj*dy
        
        if not( resname in self.res_atomlist.keys() ) :
            self.res_atomlist[resname] = []
            self.res_nomicon[resname] = set()
            self.res_nomicon[resname].add(asl)
            self.res_nomicon[resname].add(ao1)
            self.res_nomicon[resname].add(ao2)

        n = len(self.res_atomlist[resname])
        x0 = 0
        y0 = 0.5*dy
        for j in range(nj) :
            for i in range(ni) :
                n += 1
                y = y0 + j*dy
                x = x0 + i*dx + (j%2)*dx_y
                self.res_atomlist[resname].append( Atom(ao1, n, x, y, z_wall+d_so, 0.0, 0.0, 0.0) )
                self.res_atomlist[resname].append( Atom(asl, n, x, y, z_wall, 0.0, 0.0, 0.0) )
                self.res_atomlist[resname].append( Atom(ao2, n, x, y, z_wall-d_so, 0.0, 0.0, 0.0) )
                self.n_atoms += 3

    """
        Creates a flat silica monolayer grafted with alkane chains (default: hexane)
    """
    def silica_monolayer_grafted(self, z_wall, dens, n_chain=6, sp=0.450, resname1='SUB', resname2='GRA', cut_box=False) :
        
        dx = sp
        dy = sp*alpha_1
        dx_y = dx/2.0
        
        ni = int(self.box_xx/dx)
        nj = int(self.box_yy/dy)

        if cut_box :
            self.box_xx = ni*dx
            # Make sure that when cutting long y the number of layes is even (to avoid period images inconst.)
            nj += (nj%2)
            self.box_yy = nj*dy

        # Number of grafted and non-grafted silica given the surface density and the surface area [nm]
        ng = int(dens*self.box_yy*self.box_xx)
        ns = ni*nj-ng

        assert ns>=0, "Not enough quadrupoles to obtain the desired grafted sites density"
        print('Grafting '+str(ng)+' quadrupoles out of '+str(ns+ng))

        if ns>0 and not( resname1 in self.res_atomlist.keys() ) :
            self.res_atomlist[resname1] = []
            self.res_nomicon[resname1] = set()
            self.res_nomicon[resname1].add(asl)
            self.res_nomicon[resname1].add(ao1)
            self.res_nomicon[resname1].add(ao2)

        if ng>0 and not( resname2 in self.res_atomlist.keys() ) :
            self.res_atomlist[resname2] = []
            self.res_nomicon[resname2] = set()
            self.res_nomicon[resname2].add(asl)
            self.res_nomicon[resname2].add(ac1)
            self.res_nomicon[resname2].add(ao2)
            self.res_nomicon[resname2].add('C')
            self.res_nomicon[resname2].add('H1')
            self.res_nomicon[resname2].add('H2')
            self.res_nomicon[resname2].add('H3')

        n = 0
        n1 = len(self.res_atomlist[resname1])
        n2 = len(self.res_atomlist[resname2])

        idx_g = rng.choice(ng+ns,ng,replace=False)
        print(idx_g)

        x0 = 0
        y0 = 0.5*dy
        for j in range(nj) :
            for i in range(ni) :
                y = y0 + j*dy
                x = x0 + i*dx + (j%2)*dx_y
                if n in idx_g :
                    n2 += 1
                    zcs = z_wall+d_so            
                    for ic in range(n_chain-2) :
                        self.res_atomlist[resname2].append(
                            Atom('C',n2,x,y,zcs+(ic+1)*d_cc,0.0,0.0,0.0))
                        self.res_atomlist[resname2].append(
                            Atom('H1',n2,x+d_ch,y,zcs+(ic+1)*d_cc,0.0,0.0,0.0))
                        self.res_atomlist[resname2].append(
                            Atom('H2',n2,x-d_ch,y,zcs+(ic+1)*d_cc,0.0,0.0,0.0))
                    self.res_atomlist[resname2].append(
                        Atom('C',n2,x,y,zcs+(n_chain-1)*d_cc,0.0,0.0,0.0))
                    self.res_atomlist[resname2].append( 
                        Atom('H1',n2,x+d_ch,y,zcs+(n_chain-1)*d_cc,0.0,0.0,0.0))
                    self.res_atomlist[resname2].append( 
                        Atom('H2',n2,x-d_ch,y,zcs+(n_chain-1)*d_cc,0.0,0.0,0.0))
                    self.res_atomlist[resname2].append( 
                        Atom('H3',n2,x,y,zcs+(n_chain-1)*d_cc+d_ch,0.0,0.0,0.0))
                    self.res_atomlist[resname2].append( Atom(ao2, n2, x, y, z_wall-d_so, 0.0, 0.0, 0.0) )
                    self.res_atomlist[resname2].append( Atom(asl, n2, x, y, z_wall, 0.0, 0.0, 0.0) )
                    self.res_atomlist[resname2].append( Atom(ac1, n2, x, y, zcs, 0.0, 0.0, 0.0) ) 
                    self.n_atoms += (3+3*n_chain-2)
                else :
                    n1 += 1
                    self.res_atomlist[resname1].append( Atom(ao1, n1, x, y, z_wall+d_so, 0.0, 0.0, 0.0) )
                    self.res_atomlist[resname1].append( Atom(asl, n1, x, y, z_wall, 0.0, 0.0, 0.0) )
                    self.res_atomlist[resname1].append( Atom(ao2, n1, x, y, z_wall-d_so, 0.0, 0.0, 0.0) )
                    self.n_atoms += 3
                n += 1

    """
        Creates a corrugated silica monolayes at the prescribed z coordinate
    """
    def silica_monolayer_rough ( 
        self, 
        z_wall, 
        amplitude, 
        wave_number, 
        wave_offset,
        mode='lenght',
        ni=None,
        sp=0.450, 
        resname='SUB', 
        cut_box=False,
        ) :

        # Local differential geometry
        x = 0.0
        a = amplitude*wave_number
        a2 = (amplitude*wave_number)**2
        norm = lambda x : np.sqrt( 1.0 + ( a2*(np.cos(wave_number*x+wave_offset))**2 ) )
        scale_x = lambda x : 1.0 / np.sqrt( 1.0 + a2*(np.cos(wave_number*x+wave_offset))**2 )
        dx_so = lambda x : - ( d_so/norm(x) ) * (a*np.cos(wave_number*x+wave_offset))
        dz_so = lambda x : ( d_so/norm(x) )
        # Calculate spacing
        dx = sp
        dy = sp*alpha_1
        dx_y = dx/2.0
        
        if mode == 'lenght' :
            ni = int(r_rough(a2)*self.box_xx/dx)
        else :
            assert mode=='number' and not(ni==None), "Unrecognized mode or undefined number of lattice points"
        nj = int(self.box_yy/dy)

        if cut_box :
            self.box_xx = ni*dx/r_rough(a2)
            # Make sure that when cutting long y the number of layes is even (to avoid period images inconst.)
            nj += (nj%2)
            self.box_yy = nj*dy

        if not( resname in self.res_atomlist.keys() ) :
            self.res_atomlist[resname] = []
            self.res_nomicon[resname] = set()
            self.res_nomicon[resname].add(asl)
            self.res_nomicon[resname].add(ao1)
            self.res_nomicon[resname].add(ao2)

        n = len(self.res_atomlist[resname])
        x0 = 0.0
        y0 = 0.5*dy
        for j in range(nj) :
            x = x0
            # !!!
            x += (j%2)*dx_y*scale_x(x)
            for i in range(ni) :
                n += 1
                y = y0 + j*dy
                z_pert = amplitude*np.sin( wave_number*x + wave_offset ) + z_wall
                # Order is important! O1, SI, O2 
                self.res_atomlist[resname].append( Atom(ao1, n, x+dx_so(x), y, z_pert+dz_so(x), 0.0, 0.0, 0.0) )
                self.res_atomlist[resname].append( Atom(asl, n, x, y, z_pert, 0.0, 0.0, 0.0) )
                self.res_atomlist[resname].append( Atom(ao2, n, x-dx_so(x), y, z_pert-dz_so(x), 0.0, 0.0, 0.0) )
                self.n_atoms += 3
                for ii in range(inner_steps) :
                    x = x + (1.0/inner_steps)*scale_x(x)*dx

    """
        Creates a flat silica monolayer at the prescribed z coordinate
    """
    def silica_patches(self, z_wall, no_patches, mask_array, sp=0.450, res_pattern='WL', cut_box=False) :
        
        dx = sp
        dy = sp*alpha_1
        dx_y = dx/2.0
        
        ni = int(self.box_xx/dx)
        nj = int(self.box_yy/dy)

        if cut_box :
            self.box_xx = ni*dx
            # Make sure that when cutting long y the number of layes is even (to avoid period images inconst.)
            nj += (nj%2)
            self.box_yy = nj*dy
        
        resname = []
        for r in range(no_patches) :
            resname.append(res_pattern+str(r))
            if not( resname[r] in self.res_atomlist.keys() ) :
                self.res_atomlist[resname[r]] = []
                self.res_nomicon[resname[r]] = set()
                self.res_nomicon[resname[r]].add(asl)
                self.res_nomicon[resname[r]].add(ao1)
                self.res_nomicon[resname[r]].add(ao2)

        n = 0
        x0 = 0
        y0 = 0.5*dy
        for j in range(nj) :
            for i in range(ni) :
                n += 1
                y = y0 + j*dy
                x = x0 + i*dx + (j%2)*dx_y
                res = resname[mask_array[j + i*nj]]
                self.res_atomlist[res].append( Atom(ao1, n, x, y, z_wall+d_so, 0.0, 0.0, 0.0) )
                self.res_atomlist[res].append( Atom(asl, n, x, y, z_wall, 0.0, 0.0, 0.0) )
                self.res_atomlist[res].append( Atom(ao2, n, x, y, z_wall-d_so, 0.0, 0.0, 0.0) )
                self.n_atoms += 3
    
    def silica_checker(self, z_wall, width_x, width_y, off_x=0, off_y=0, sp=0.450, res_pattern='WL', cut_box=False) :

        dx = sp
        dy = sp*alpha_1
        dx_y = dx/2.0

        ni = int(self.box_xx/dx)
        nj = int(self.box_yy/dy)

        if cut_box :
            self.box_xx = ni*dx
            # Make sure that when cutting long y the number of layes is even (to avoid period images inconst.)
            nj += (nj%2)
            self.box_yy = nj*dy

        resname = []
        for r in range(2) :
            resname.append(res_pattern+str(r+1))
            if not( resname[r] in self.res_atomlist.keys() ) :
                self.res_atomlist[resname[r]] = []
                self.res_nomicon[resname[r]] = set()
                self.res_nomicon[resname[r]].add(asl)
                self.res_nomicon[resname[r]].add(ao1)
                self.res_nomicon[resname[r]].add(ao2)

        n = 0
        x0 = 0
        y0 = 0.5*dy
        for j in range(nj) :
            for i in range(ni) :
                n += 1
                y = y0 + j*dy
                x = x0 + i*dx + (j%2)*dx_y
                if int((x+off_x)/width_x)%2 == 0 and int((y+off_y)/width_y)%2 == 0 :
                    res = resname[0]
                else :
                    res = resname[1]
                self.res_atomlist[res].append( Atom(ao1, n, x, y, z_wall+d_so, 0.0, 0.0, 0.0) )
                self.res_atomlist[res].append( Atom(asl, n, x, y, z_wall, 0.0, 0.0, 0.0) )
                self.res_atomlist[res].append( Atom(ao2, n, x, y, z_wall-d_so, 0.0, 0.0, 0.0) )
                self.n_atoms += 3


    def silica_circle(self, xc, yc, R, new_resname='WLN', old_resname='SUB') :

        self.res_atomlist[new_resname] = []
        self.res_nomicon[new_resname] = self.res_nomicon[old_resname]
        
        n_atoms_mol = len(self.res_nomicon[old_resname])
        auxillary_atom_list = []
        auxillary_molecule = []

        n = 0
        r = 0
        for a in self.res_atomlist[old_resname] :
            r = max(r, np.sqrt( (a.pos_x-xc)**2 + (a.pos_y-yc)**2 ) )
            auxillary_molecule.append(a)
            n += 1
            if n == n_atoms_mol :
                if r <= R :
                    self.res_atomlist[new_resname].extend(auxillary_molecule)
                else :
                    auxillary_atom_list.extend(auxillary_molecule)
                auxillary_molecule = []
                n = 0
                r = 0

        self.res_atomlist[old_resname] = auxillary_atom_list


    def silica_cristobalite(self) :
        pass    


    def silica_amorphous_read(self, file_name) :

        n_lines = count_line(file_name)
        n = 0 
        in_file = open(file_name, 'r')
        self.header = "Amorphous silica slab"
        self.res_atomlist['SLB'] = []
        self.res_nomicon['SLB'] = set()
        for line in in_file :
            n += 1
            if n == n_lines :
                cols = line.split()
                self.box_xx = float(cols[0])
                self.box_yy = float(cols[1])
                self.box_zz = float(cols[2])
            elif n == 2 :
                self.n_atoms = int(line)
            elif n > 2 :
                line_data = read_gro_line(line, read_vel=False)
                if line_data[2] == '1' :
                    atom_name = 'Si'
                else :
                    atom_name = 'O'
                self.res_atomlist['SLB'].append( Atom(atom_name, line_data[3], \
                    line_data[4], line_data[5], line_data[6], 0.0, 0.0, 0.0) )
                self.res_nomicon['SLB'].add(atom_name)
        in_file.close()

    """
        To produce bulk silica topologies (to check rdf, for example)
    """
    def silica_amorphous_bulk(self) : 

        h_list = list(self.res_atomlist.keys())

        self.res_atomlist['SBS'] = []
        self.res_nomicon['SBS'] = {'BS'}

        self.res_atomlist['SBO'] = []
        self.res_nomicon['SBO'] = {'BO'}

        for h in h_list :
            for a in self.res_atomlist[h] :

                if a.label == "Si" :
                    a.label='BS'
                    self.res_atomlist['SBS'].append(a)
                elif a.label == "O" :
                    a.label='BO'
                    self.res_atomlist['SBO'].append(a)

            self.res_atomlist.pop(h)
            self.res_nomicon.pop(h)

    """
        To produce a surface
    """
    def silica_amorphous_cut(self, z_inf, z_half, z_sup, resname='SLB') :

        n_surf = 0
        surf_atoms = []
        n_subs = 0
        subs_atoms = []

        for a in self.res_atomlist[resname] :
            if a.pos_z > z_inf and a.pos_z <= z_half :
                subs_atoms.append(a)
                n_subs += 1
            if a.pos_z > z_half and a.pos_z <= z_sup :
                surf_atoms.append(a)
                n_surf += 1

        self.res_atomlist.pop(resname)
        self.res_nomicon.pop(resname)

        self.res_atomlist['SUR'] = surf_atoms
        self.res_nomicon['SUR'] = {'O', 'Si', 'DO', 'DSi'}
        self.res_atomlist['SUB'] = subs_atoms
        self.res_nomicon['SUB'] = {'O', 'Si'}
        self.n_atoms = n_subs+n_surf

    def silica_amorphous_carve_sine(self, z0, h, w, phi=0, cut_sub=False) :

        n_surf = len(self.res_atomlist['SUR'])
        n_surf_new = 0
        surf_atoms_new = []

        z_sine = lambda x : z0 + h*np.sin(w*x+phi)

        for a in self.res_atomlist['SUR'] :
            z = a.pos_z
            x = a.pos_x
            if z <= z_sine(x) :
                surf_atoms_new.append(a)
                n_surf_new += 1
                
        self.res_atomlist['SUR'] = surf_atoms_new
        self.n_atoms -= (n_surf-n_surf_new)

        if cut_sub :
            n_sub = len(self.res_atomlist['SUB'])
            n_sub_new = 0
            sub_atoms_new = []
            for a in self.res_atomlist['SUB'] :
                z = a.pos_z
                x = a.pos_x
                if z > z_sine(x) -z0 + h :
                    sub_atoms_new.append(a)
                    n_sub_new += 1
            self.res_atomlist['SUB'] = sub_atoms_new
            self.n_atoms -= (n_sub-n_sub_new)


    def silica_amorphous_carve_pillar(self, z0, h, b, s=0) :

        n_surf = len(self.res_atomlist['SUR'])
        n_surf_new = 0
        surf_atoms_new = []

        z_pillar = lambda x : (z0+h)*(int((x+s)/b)%2==0) + (z0-h)*(int((x+s)/b)%2==1)

        for a in self.res_atomlist['SUR'] :
            z = a.pos_z
            x = a.pos_x
            if z <= z_pillar(x) :
                surf_atoms_new.append(a)
                n_surf_new += 1

        self.res_atomlist['SUR'] = surf_atoms_new
        self.n_atoms -= (n_surf-n_surf_new)

    def silica_amorphous_carve_triangle(self, z0, h, c, b, s=0) :

        n_surf = len(self.res_atomlist['SUR'])
        n_surf_new = 0
        surf_atoms_new = []

        xi = lambda x : x-int(x/(c+b))*(c+b)
        z_triangle = lambda sx : z0 + (sx<c)*(h/c)*sx + (sx>=c)*(h-(h/b)*(sx-c))

        for a in self.res_atomlist['SUR'] :
            z = a.pos_z
            x = a.pos_x
            sx = xi(x)
            if z <= z_triangle(sx) :
                surf_atoms_new.append(a)
                n_surf_new += 1

        self.res_atomlist['SUR'] = surf_atoms_new
        self.n_atoms -= (n_surf-n_surf_new)


    def silica_amorphous_carve_pillar_3D(self,hx,hy,hz,lx,ly) :

        n_surf = len(self.res_atomlist['SUR'])
        n_surf_new = 0
        surf_atoms_new = []

        dx = hx+lx
        dy = hy+ly

        for a in self.res_atomlist['SUR'] :
            x = a.pos_x
            y = a.pos_y
            z = a.pos_z
            if z <= self.box_zz-hz :
                surf_atoms_new.append(a)
                n_surf_new += 1
            else :
                xloc = x-dx*int(x/dx)
                yloc = y-dy*int(y/dy)
                print('xloc ',xloc)
                print('yloc', yloc)
                if xloc>=lx and yloc>=ly :
                    surf_atoms_new.append(a)
                    n_surf_new += 1


        self.res_atomlist['SUR'] = surf_atoms_new
        self.n_atoms -= (n_surf-n_surf_new)


    """
    Remove the excess Si or O  strating from the atoms having the smaller z coordinate
    """
    def silica_amorphous_equalize_from_below(self) :

        ZLARGE = 10*self.box_zz

        nsi = sum(a.label == "Si" for a in self.res_atomlist['SUB'])
        nsi += sum(a.label == "Si" for a in self.res_atomlist['SUR'])
        no = sum(a.label == "O" for a in self.res_atomlist['SUB'])
        no += sum(a.label == "O" for a in self.res_atomlist['SUR'])
        balance = 2*nsi-no
        print("balance before = "+str(balance))

        if balance > 0 :
            # Too many silicon atoms
            if balance%2>0 :
                # Remove a oxygen atom
                amin = min(self.res_atomlist['SUB'], key = lambda a: a.pos_z+ZLARGE*(a.label=='Si'))
                self.res_atomlist['SUB'].remove(amin)
                self.n_atoms -= 1
                balance += 1
            # Remove "balance/2" silicon atoms
            while balance>0 :
                amin = min(self.res_atomlist['SUB'], key = lambda a: a.pos_z+ZLARGE*(a.label=='O'))
                self.res_atomlist['SUB'].remove(amin)
                self.n_atoms -= 1
                balance -= 2
        elif balance < 0 :
            # Too many oxygen atoms, remove "balance" oxgen atoms
            while balance<0 :
                amin = min(self.res_atomlist['SUB'], key = lambda a: a.pos_z+ZLARGE*(a.label=='Si'))
                self.res_atomlist['SUB'].remove(amin)
                self.n_atoms -= 1
                balance += 1

        # nsi = sum(a.label == "Si" for a in self.res_atomlist['SUB'])
        # nsi += sum(a.label == "Si" for a in self.res_atomlist['SUR'])
        # no = sum(a.label == "O" for a in self.res_atomlist['SUB'])
        # no += sum(a.label == "O" for a in self.res_atomlist['SUR'])
        # balance = 2*nsi-no
        print("balance now = "+str(balance))

    """
    Converg back to the atomtypes of LAMMPS tutorial
    """
    def silica_amorphous_gro2lam(self) :

        sil_atoms = []
        list_of_keys = list(self.res_atomlist.keys())
        print(list_of_keys)

        for key in list_of_keys :
            for a in self.res_atomlist[key] :
                if a.label == 'Si' or a.label == 'DSi' :
                    a.label = 1
                    sil_atoms.append(a)
                if a.label == 'O' or a.label == 'DO' :
                    a.label = 2
                    sil_atoms.append(a)
            self.res_atomlist.pop(key)
            self.res_nomicon.pop(key)

        self.res_atomlist['SIL'] = sil_atoms
        self.res_nomicon['SIL'] = {'1', '2'}

        list_of_keys = list(self.res_atomlist.keys())
        print(list_of_keys)

    """
    Go back to the notes about structure functions! 
    def silica_amorphous_make_rough(self, ...)
    !!! THIS FUNCTION IS INCORRECT AS IT DOES NOT ACCOUNT FOR PBC !!!
    """
    def silica_amorphous_connect(self, r_bond=0.2) :

        self.res_nomicon['SUR'].add('DO')
        self.res_nomicon['SUR'].add('DSi')

        n_surf = len(self.res_atomlist['SUR'])
        n = 0

        for n in range(n_surf) :
            a = self.res_atomlist['SUR'][n]
            print("Atom "+str(n+1)+"/"+str(n_surf))
            if a.label == 'O' :
                nt = 2
                tag = 'DO'
                tar = 'Si'
            elif a.label == 'Si' :
                nt = 4
                tag = 'DSi'
                tar = 'O'
            n_neigh = 0
            # Loop over SUR
            for b in self.res_atomlist['SUR'] :
                if n_neigh >= nt :
                    break
                elif b.label == tar :
                    d = self.pcb_dist(b,a)
                    if d <= r_bond :
                        n_neigh += 1
            # Loop over SUB
            for b in self.res_atomlist['SUB'] :
                if n_neigh >= nt :
                    break
                elif b.label == tar :
                    d = self.pcb_dist(b,a)
                    if d <= r_bond :
                        n_neigh += 1
            if n_neigh < nt :
                print("Found a dandling atom!")
                self.res_atomlist['SUR'][n].label = tag

    """
        Relabel atoms and residues so that every specie is a different molecule
    """
    def silica_amorphous_relabel(self) :

        self.res_atomlist['SRS'] = []
        self.res_nomicon['SRS'] = {'RS'}

        self.res_atomlist['SRO'] = []
        self.res_nomicon['SRO'] = {'RO'}

        self.res_atomlist['DAS'] = []
        self.res_nomicon['DAS'] = {'DS'}

        self.res_atomlist['DAO'] = []
        self.res_nomicon['DAO'] = {'DO'}

        self.res_atomlist['SBS'] = []
        self.res_nomicon['SBS'] = {'BS'}

        self.res_atomlist['SBO'] = []
        self.res_nomicon['SBO'] = {'BO'}

        for a in self.res_atomlist['SUR'] :
            if a.label == "Si" :
                a.label='RS'
                self.res_atomlist['SRS'].append(a)
            elif a.label == "O" :
                a.label='RO'
                self.res_atomlist['SRO'].append(a)
            elif a.label == "DSi" :
                a.label='DS'
                self.res_atomlist['DAS'].append(a)
            else :
                a.label='DO'
                self.res_atomlist['DAO'].append(a)
        
        for a in self.res_atomlist['SUB'] :
            if a.label == "Si" :
                a.label='BS'
                self.res_atomlist['SBS'].append(a)
            elif a.label == "O" :
                a.label='BO'
                self.res_atomlist['SBO'].append(a)

        self.res_atomlist.pop('SUB')
        self.res_atomlist.pop('SUR')
        self.res_nomicon.pop('SUB')
        self.res_nomicon.pop('SUR')


    """
        BETA: add a sinal group for each dandling silicon atom
    """
    def silica_add_silanol(self) :

        # Number of silanol groups to be added
        ns = len(self.res_atomlist['DAS'])
        
        # Uniform random angle
        random_phi = 2.0*np.pi*rng.random(ns)
        scaling_factor = 1.0 / np.sqrt( 1.0 + np.tan(random_phi)*np.tan(random_phi) )
        sign = -2 * (random_phi>0.5*np.pi) * (random_phi<=1.5*np.pi) + 1

        # Projected distance of hydrogen atoms
        """
        Names:
        d_ct2_oh1 = 0.142
        d_oh1_h = 0.096
        a_ct2_oh1_h = 106.0000
        """
        dzh = d_oh1_h*np.sin( np.deg2rad(a_ct2_oh1_h-90) )
        distance_factor = d_oh1_h*np.cos( np.deg2rad(a_ct2_oh1_h-90) )
        dxh = sign*distance_factor*scaling_factor
        dyh = dxh*np.tan(random_phi)

        # Create new residue SIL
        self.res_atomlist['SIL'] = []
        self.res_nomicon['SIL'] = {'Si', 'O', 'H'}

        nsil = self.n_atoms
        i = 0
        for a in self.res_atomlist['DAS'] :
            nsil += 1
            a.label = 'Si'
            a.n = nsil
            self.res_atomlist['SIL'].append(a)
            nsil += 1
            self.res_atomlist['SIL'].append(Atom( 'O', nsil, a.pos_x, a.pos_y, a.pos_z+d_ct2_oh1 ))
            nsil += 1
            self.res_atomlist['SIL'].append(Atom( 'H', nsil, a.pos_x+dxh[i], a.pos_y+dyh[i], a.pos_z+dzh+d_ct2_oh1 ))
            i += 1

        self.res_atomlist.pop('DAS')
        self.res_nomicon.pop('DAS')

        self.n_atoms += 2*ns


    """
        Carves a sphere of molecules of a given residue
    """
    def carve_sphere(self, radius, cx=None, cy=None, cz=None, resname='SOL', n_atom_mol=None) :
       
        print("Carving a 3D sphere")
        
        cut2 = radius**2
        if cx == None :
            cx = 0.5*self.box_xx
        if cy == None :
            cy = 0.5*self.box_yy
        if cz == None :
            cz = 0.5*self.box_zz

        if not( resname in self.res_nomicon.keys() ) :
            print("Specified residue not present in current configuration")

        else :
            new_atom_list = []
            if n_atom_mol == None :
                n_atom_mol = len(self.res_nomicon[resname])
            loc_atom_list = [None] * n_atom_mol
            n_count = 0
            n = 0
            for atom in self.res_atomlist[resname] :
                x = atom.pos_x
                y = atom.pos_y
                z = atom.pos_z
                loc_atom_list[n_count] = atom
                if n_count == 0 :
                    r2 = (x-cx)**2 + (y-cy)**2 + (z-cz)**2
                else :
                    r2 = max(r2, (x-cx)**2 + (y-cy)**2 + (z-cz)**2)
                if n_count == (n_atom_mol-1) and r2 < cut2 :
                    for k in range(n_atom_mol) :
                        new_atom_list.append(loc_atom_list[k])
                        n = n + 1
                n_count = ( n_count + 1 ) % n_atom_mol
            self.n_atoms -= (len(self.res_atomlist[resname])-n)
            self.res_atomlist[resname] = new_atom_list
    
    """
        Carves a 3D cylinder (quasi-2d droplets)
    """
    def carve_cylinder(self, radius, cx=None, cz=None, resname='SOL', n_atom_mol=None) :
        
        print("Carving a cylinder (2D droplet)")
        
        cut2 = radius**2
        if cx == None :
            cx = 0.5*self.box_xx
        if cz == None :
            cz = 0.5*self.box_zz

        if not( resname in self.res_nomicon.keys() ) :
            print("Specified residue not present in current configuration")

        else :
            new_atom_list = []
            if n_atom_mol == None :
                n_atom_mol = len(self.res_nomicon[resname])
            loc_atom_list = [None] * n_atom_mol
            n_count = 0
            n = 0
            for atom in self.res_atomlist[resname] :
                x = atom.pos_x
                z = atom.pos_z
                loc_atom_list[n_count] = atom
                if n_count == 0 :
                    r2 = (x-cx)**2 + (z-cz)**2
                else :
                    r2 = max(r2, (x-cx)**2 + (z-cz)**2)
                if n_count == (n_atom_mol-1) and r2 < cut2 :
                    for k in range(n_atom_mol) :
                        new_atom_list.append(loc_atom_list[k])
                        n = n + 1
                n_count = ( n_count + 1 ) % n_atom_mol
            self.n_atoms -= (len(self.res_atomlist[resname])-n)
            self.res_atomlist[resname] = new_atom_list

    """
        Carves a rectangle (to inizialize confined simulations)
    """
    def carve_rectangle(self, crop_x, crop_z, crop_y=None, resname='SOL', n_atom_mol=None) :

        if not( resname in self.res_nomicon.keys() ) :
            print("Specified residue not present in current configuration")

        else :
            new_atom_list = []
            if n_atom_mol == None :
                n_atom_mol = len(self.res_nomicon[resname])
            loc_atom_list = [None] * n_atom_mol
            n_count = 0
            n = 0

            if crop_y == None :

                for atom in self.res_atomlist[resname] :
                    x = atom.pos_x
                    y = atom.pos_y
                    z = atom.pos_z
                    loc_atom_list[n_count] = atom
                    if n_count == 0 :
                        in_rect = (x > crop_x[0]) * ( x <= crop_x[1] ) \
                            * (z > crop_z[0]) * ( z <= crop_z[1] )
                    else :
                        in_rect = in_rect * (x > crop_x[0]) * ( x <= crop_x[1] ) \
                            * (z > crop_z[0]) * ( z <= crop_z[1] )
                    if n_count == (n_atom_mol-1) and in_rect :
                        for k in range(n_atom_mol) :
                            new_atom_list.append(loc_atom_list[k])
                            n = n + 1
                    n_count = ( n_count + 1 ) % n_atom_mol
                self.n_atoms -= (len(self.res_atomlist[resname])-n)
                self.res_atomlist[resname] = new_atom_list

            else :

                for atom in self.res_atomlist[resname] :
                    x = atom.pos_x
                    y = atom.pos_y
                    z = atom.pos_z
                    loc_atom_list[n_count] = atom
                    if n_count == 0 :
                        in_rect = (x > crop_x[0]) * ( x <= crop_x[1] ) \
                            * (z > crop_z[0]) * ( z <= crop_z[1] ) \
                            * (y > crop_y[0]) * ( y <= crop_y[1] )
                    else :
                        in_rect = in_rect * (x > crop_x[0]) * ( x <= crop_x[1] ) \
                            * (z > crop_z[0]) * ( z <= crop_z[1] ) \
                            * (y > crop_y[0]) * ( y <= crop_y[1] )
                    if n_count == (n_atom_mol-1) and in_rect :
                        for k in range(n_atom_mol) :
                            new_atom_list.append(loc_atom_list[k])
                            n = n + 1
                    n_count = ( n_count + 1 ) % n_atom_mol
                self.n_atoms -= (len(self.res_atomlist[resname])-n)
                self.res_atomlist[resname] = new_atom_list

    """
        Resize to adapt to the period box
    """
    def resize_to_box(self, resname='SOL', n_atom_mol=None) :
        
        if not( resname in self.res_nomicon.keys() ) :
            print("Specified residue not present in current configuration")

        else :
            new_atom_list = []
            if n_atom_mol == None :
                n_atom_mol = len(self.res_nomicon[resname])
            loc_atom_list = [None] * n_atom_mol
            n_count = 0
            n = 0
            for atom in self.res_atomlist[resname] :
                x = atom.pos_x
                y = atom.pos_y
                z = atom.pos_z
                loc_atom_list[n_count] = atom
                if n_count == 0 :
                    in_box = (x > 0) * ( x <= self.box_xx ) \
                        * (y > 0) * ( y <= self.box_yy ) \
                        * (z > 0) * ( z <= self.box_zz )
                else :
                    in_box = in_box * (x > 0) * ( x <= self.box_xx ) \
                        * (y > 0) * ( y <= self.box_yy ) \
                        * (z > 0) * ( z <= self.box_zz )
                if n_count == (n_atom_mol-1) and in_box :
                    for k in range(n_atom_mol) :
                        new_atom_list.append(loc_atom_list[k])
                        n = n + 1
                n_count = ( n_count + 1 ) % n_atom_mol
            self.n_atoms -= (len(self.res_atomlist[resname])-n)
            self.res_atomlist[resname] = new_atom_list

    """
        Merges with another configuration
    """
    def merge(self, conf) :
        
        self.n_atoms += conf.n_atoms
        
        offset_x = 0.5 * (self.box_xx-conf.box_xx)
        offset_y = 0.5 * (self.box_yy-conf.box_yy)
        offset_z = 0.5 * (self.box_zz-conf.box_zz)

        for resname in conf.res_nomicon.keys() :
            if not( resname in self.res_nomicon.keys() ) :
                self.res_atomlist[resname] = []
                self.res_nomicon[resname] = set()
                for atom in conf.res_atomlist[resname] :
                    self.res_atomlist[resname].append( atom )
                    self.res_nomicon[resname].add( atom.label )
            else :
                n_off = len(self.res_atomlist[resname])
                for atom in conf.res_atomlist[resname] :
                    atom.index += n_off
                    atom.pos_x += offset_x
                    atom.pos_y += offset_y
                    atom.pos_z += offset_z
                    self.res_atomlist[resname].append( atom )
    
    """
        Shifts residues
    """
    def shift(self, dx, dy, dz, resname='SOL') :
        
        for k in range(len(self.res_atomlist[resname])) :
            self.res_atomlist[resname][k].pos_x += dx
            self.res_atomlist[resname][k].pos_y += dy
            self.res_atomlist[resname][k].pos_z += dz

    """
        Shift the same residue, but in opposite direction depending whether is on top or on bottom
    """
    def shift_walls(self, dx, dy, dz, resname='SUB') :
        
        for k in range(len(self.res_atomlist[resname])) :
            if self.res_atomlist[resname][k].pos_z > 0.5*self.box_zz :
                self.res_atomlist[resname][k].pos_x += dx
                self.res_atomlist[resname][k].pos_y += dy
                self.res_atomlist[resname][k].pos_z += dz
            else :
                self.res_atomlist[resname][k].pos_x -= dx
                self.res_atomlist[resname][k].pos_y -= dy
                self.res_atomlist[resname][k].pos_z -= dz

    """
        Reflect and shift (used for amorphous silica substrates)
    """
    def reflect_walls(self, dx=0, resname='SUB') :
        
        assert (dx>=0), "Negative horizontal shift unsupported"

        for k in range(len(self.res_atomlist[resname])) :
            x = self.res_atomlist[resname][k].pos_x
            new = ((x+dx)/self.box_xx-int((x+dx)/self.box_xx))*self.box_xx
            self.res_atomlist[resname][k].pos_x = new
            z = self.res_atomlist[resname][k].pos_z
            self.res_atomlist[resname][k].pos_z = self.box_zz-z

    """
        Self explanatory
    """
    def remove_residue(self, resname) :

        nres = len(self.res_atomlist[resname])
        self.res_atomlist.pop(resname)
        self.res_nomicon.pop(resname)
        self.n_atoms -= nres
    
    """
        Add uniform constant velocity to a residue
    """
    def add_uniform_velocity(self, dvx, dvy, dvz, resname='SOL') :

        for k in range(len(self.res_atomlist[resname])) :
            self.res_atomlist[resname][k].vel_x += dvx
            self.res_atomlist[resname][k].vel_y += dvy
            self.res_atomlist[resname][k].vel_z += dvz

    """
        Assign (!) a linear velocity profile (with zero mean) along z in x direction to a residue
    """
    def assign_shear_velocity(self, deform_rate, resname='SOL') :
        
        bz = self.box_zz
        for k in range(len(self.res_atomlist[resname])) :
            z = self.res_atomlist[resname][k].pos_z
            if z < 0 :
                z = bz+z
            if z >= bz :
                z = z-bz
            # self.res_atomlist[resname][k].vel_x += deform_rate*(z-0.5*bz)/bz
            self.res_atomlist[resname][k].vel_x = deform_rate*(z-0.5*bz)/bz

    """
        Bin velocity profile in x direction along z
    """
    def bin_velocity_profile(self, n_bins=10) :

        bz = self.box_zz
        dz = bz/n_bins
        z_range = np.linspace(dz, bz-dz, n_bins)
        x_velocity = np.zeros(n_bins)
        n_count = np.zeros(n_bins)

        for resname in self.res_atomlist.keys() :

            for k in range(len(self.res_atomlist[resname])) :

                z = self.res_atomlist[resname][k].pos_z
                if z < 0 :
                    z = bz+z
                if z >= bz :
                    z = z-bz
                v = self.res_atomlist[resname][k].vel_x
                idx = int(z/dz)
                x_velocity[idx] += v
                n_count[idx] += 1

        x_velocity /= n_count

        return z_range, x_velocity

    """
        Print informations
    """
    def print_info(self) :
        print("### MDCONF ###")
        print(self.header)
        print("no. atoms = "+str(self.n_atoms))
        print("Lx = "+str(self.box_xx))
        print("Ly = "+str(self.box_yy))
        print("Lz = "+str(self.box_zz))
        for resname in self.res_nomicon.keys() :
            nm = int(len(self.res_atomlist[resname])/len(self.res_nomicon[resname]))
            print("Res: "+resname+", #molecules = "+str(nm)+", atom types "+str(self.res_nomicon[resname]))
        print("##############")












