import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.special

# DEFINITION THAT SHOULD NOT BE CHANGED BY THE USER #

# Silicium-oxygen distance in silica quadrupoles
d_so = 0.151
# Atom types in silica quadrupoles
asl = "SI"
ao1 = "O1"
ao2 = "O2"
# Hexagonal silica lattice parameters
alpha_1 = np.sqrt(3.0/4.0)
alpha_2 = np.sqrt(2.0/3.0)
# NUmber of steps to differnciate local geometry on rough substrates
inner_steps = 5
# Effective surface as function of roughness parameter
r_rough = lambda a2 : (2.0/np.pi) * np.sqrt(a2+1.0) * sc.special.ellipe(a2/(a2+1.0))

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

def read_gro_line( line ) :
    line_data = [None]*10
    line_data[0] = int(line[0:5].strip())                   # Residue serial
    line_data[1] = str(line[5:10].strip())                  # Residue name
    line_data[2] = str(line[10:15].strip())                 # Atom name
    line_data[3] = int(line[15:20].strip())                 # Atom serial
    line_data[4] = float(line[20:28].strip())               # Pos X
    line_data[5] = float(line[28:36].strip())               # Pos Y
    line_data[6] = float(line[36:44].strip())               # Pos Z
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
        Creates a corrigated silica monolayes at the prescribed z coordinate
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
    
    """
        Carves a sphere of molecules of a given residue
    """
    def carve_sphere(self, radius, cx=None, cy=None, cz=None, resname='SOL') :
       
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
        Carves a rectangle (to inizialize confined simulations)
    """
    def carve_rectangle(self, crop_x, crop_z, resname='SOL') :
        
        if not( resname in self.res_nomicon.keys() ) :
            print("Specified residue not present in current configuration")

        else :
            new_atom_list = []
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
        Add linear velocity to a residue
    """
    def add_velocity(self, dvx, dvy, dvz, resname='sol') :

        for k in range(len(self.res_atomlist[resname])) :
            self.res_atomlist[resname][k].vel_x += dvx
            self.res_atomlist[resname][k].vel_y += dvy
            self.res_atomlist[resname][k].vel_z += dvz

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












