from math import sqrt
from math import sin
from math import cos

import numpy as np

"""
    Library for creating and manipulating MD configurations (reference: GROMACS)
    MAKE THIS MORE OBJECT-ORIENTED!
"""

"""
    FUNCTION lj_substrate
    Generates a Lennard-Jones substrate; the arguments refer to the number of
    atoms in each direction, the spacing, the angles (?) between directions of
    the crystal structure and the output file name.
    THIS FUNCTION IS NOT COMPLETE! Need implementation for a generic kind of
    crystalline structure; at the present moment FCC (?) is default.
"""
def lj_substrate (
    ni = 10,
    nj = 10,
    nk = 10,
    sp = 2.7,
    alpha_1 = sqrt(3.0/4.0),
    alpha_2 = sqrt(2.0/3.0),
    file_name = 'lj_substrate.pdb' ):

    dx = sp
    dy = sp*alpha_1
    dz = sp*alpha_2

    lj_file = open(file_name, 'w')

    lj_file.write( "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" %
       (ni*dx,nj*dy,nk*dz,90.0,90.0,90.0) );

    # THESE NEED TO BE MODIFIED ACCORDINGLY W.R.T. ALPHA1 AND ALPHA2
    ###
    dx_y = dx/2;
    dx_z = dx/2;
    dy_z = dy/3;
    x0 = dx/4;
    y0 = dy/6;
    z0 = dz/2;
    ###

    n = 1;
    for k in range(nk):
        if k < nk-1 :
            a = "CUB"
        else:
            a = "CUS"
        z = z0 + k*dz
        for i in range(ni):
            for j in range(nj):
                y = y0 + j*dy + k*dy_z
                x = x0 + i*dx + j*dx_y + k*dx_z
                lj_file.write( "%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM",n % 100000,a,"SUB",1,x,y,z) )
                n = n+1

    lj_file.close()


"""
    FUNCTION quad_substrate
    Generates a SiO2 substrate
    [...]
"""
def quad_substrate (
    ni = 10,
    nj = 10,
    nk = 1,
    file_name = 'quad_substrate.pdb' ):

    # Lattice parameters
    sp = 4.50
    alpha_1 = sqrt(3.0/4.0)
    alpha_2 = sqrt(2.0/3.0)

    # Atom types
    asl = "SI"
    ao1 = "O1"
    ao2 = "O2"

    # Distance of oxigen atoms
    d_so = 1.51

    # Calculate spacing
    dx = sp;
    dy = sp*alpha_1;
    dz = sp*alpha_2;

    # Row spacing
    dx_y = dx/2.0
    dx_z = 0.0
    dy_z = 0.0

    # Starting positions
    x0 = dx/4.0
    y0 = dy/6.0
    z0 = 3.0

    quad_file = open(file_name, 'w')

    quad_file.write( "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" %
       (ni*dx,nj*dy,nk*dz,90.0,90.0,90.0) );

    n = 1
    # Loop over all layers
    for k in range(nk) :
    # Check if bulk or not
        if k < (nk - 1) :
            a = "CUB"
        else :
            a = "CUS"
        z = z0 + k*dz
        # Loop over rows and columns
        for i in range(ni) :
            for j in range(nj) :
                y = y0 + j*dy + k*dy_z
                x = x0 + i*dx + j*dx_y + k*dx_z
                quad_file.write("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM", n % 100000, ao1, "SUB", 1, x, y, z+d_so) )
                quad_file.write("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM", n % 100000, asl, "SUB", 1, x, y, z) )
                quad_file.write("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM", n % 100000, ao2, "SUB", 1, x, y, z-d_so) )
                n = n+1

    quad_file.close()


"""
    FUNCTION quad_substrate_wave
    Generates a sinusoidal-shaped SiO2 substrate
    [...]
"""
def quad_substrate_wave (
    amplitude,
    amplitude_offset,
    wave_number,
    wave_offset,
    bend,
    ni = 10,
    nj = 10,
    nk = 1,
    file_name = 'quad_substrate_wave.pdb'
    ):

    # Lattice parameters
    # sp = 4.50
    sp = 4.50/sqrt(2)
    alpha_1 = sqrt(3.0/4.0)
    alpha_2 = sqrt(2.0/3.0)

    # Atom types
    asl = "SI"
    ao1 = "O1"
    ao2 = "O2"

    x = 0.0
    norm = lambda : sqrt( 1.0 + 1.0/(amplitude**2*wave_number**2*cos(wave_number*x+wave_offset)**2) )

    d_so = 1.51

    # Distance of oxigen atoms
    if bend == False :
        dx_so = lambda : 0.00
        dz_so = lambda : d_so
    else :
        dx_so = lambda : ( d_so/norm() )
        dz_so = lambda : - ( d_so/norm() ) / (amplitude*wave_number*cos(wave_number*x+wave_offset))

    # Calculate spacing
    dx = sp;
    dy = sp*alpha_1;
    dz = sp*alpha_2;

    # Row spacing
    dx_y = dx/2.0
    dx_z = 0.0
    dy_z = 0.0

    # Starting positions
    x0 = dx/4.0
    y0 = dy/6.0
    z0 = 3.0

    quad_wave_file = open(file_name, 'w')

    quad_wave_file.write( "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" %
       (ni*dx,nj*dy,nk*dz,90.0,90.0,90.0) );

    n = 1
    # Loop over all layers
    for k in range(nk) :
    # Check if bulk or not
        if k < (nk - 1) :
            a = "CUB"
        else :
            a = "CUS"
        z = z0 + k*dz
        # Loop over rows and columns
        for i in range(ni) :
            for j in range(nj) :
                y = y0 + j*dy + k*dy_z
                x = x0 + i*dx + j*dx_y + k*dx_z
                z_pert = z + amplitude*sin( wave_number*x + wave_offset ) + amplitude_offset
                assert z_pert>=0, "Negative height"
                quad_wave_file.write("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM", n % 100000, ao1, "SUB", 1, x+dx_so(), y, z_pert+dz_so()) )
                quad_wave_file.write("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM", n % 100000, asl, "SUB", 1, x, y, z_pert) )
                quad_wave_file.write("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                    ("ATOM", n % 100000, ao2, "SUB", 1, x-dx_so(), y, z_pert-dz_so()) )
                n = n+1

    quad_wave_file.close()


"""
    FUNCTION carve_2D_droplet
    Carves a 2D cylindric droplet out of a sufficienly large fluid cell; the
    function works for water, maybe in future it will see a generalization to
    an arbitrary fluid.
"""
def carve_2D_droplet(
    carving_factor,
    input_file_name,
    output_file_name = 'droplet.pdb',
    mode = 'p'
    ):

    cut2 = 0.0
    cx = 0.0
    cz = 0.0

    # Dummy
    n = 0

    # Water molecule
    n_atom_mol = 3
    n_count = 0

    liquid_file = open(input_file_name, 'r')
    droplet_file = open(output_file_name, 'w+')

    with liquid_file as f :
        for line in f :
            cols = line.split()
            if mode == 'p' and cols[0] == "CRYST1":
                cx = float(cols[1])/2.0
                cz = float(cols[3])/2.0
                cut2 = carving_factor**2*cz**2
                print("r = %f, cx = %f, cz = %f" % (sqrt(cut2), cx, cz) )
                droplet_file.write(line)
                break
        for line in f :
            cols = line.split()
            if cols[0] == "ATOM" :
                n_count = ( n_count + 1 ) % n_atom_mol
                if n_count == 1 :
                    x = float(cols[5])
                    z = float(cols[7])
                    r2 = (x-cx)**2 + (z-cz)**2
                if r2 < cut2 :
                    droplet_file.write(line)
                    n = n + 1
        print("n atoms = %d" % n)

    droplet_file.close()
    liquid_file.close()


"""
    FUNCTION adapt_to_2D_subtrate
    [...]
"""
def merge_to_substrate(
    substrate_file,
    liquid_file,
    system_file = 'system.pdb'
    ):

    sub = open(substrate_file, 'r')
    liq = open(liquid_file, 'r')
    sys = open(system_file, 'w+')

    for line in sub :
        cols = line.split()
        if cols[0] == "CRYST1" :
            crystal_col = cols
            break

    sub_max_y = float(crystal_col[2])

    for line in liq :
        cols = line.split()
        if cols[0] == "CRYST1" :
            crystal_col[3] = cols[3]
            break

    print(crystal_col)

    # Writing box parameters
    sys.write((6-len(crystal_col[0]))*' '+crystal_col[0])
    sys.write((15-6-len(crystal_col[1]))*' '+crystal_col[1])
    sys.write((24-15-len(crystal_col[2]))*' '+crystal_col[2])
    sys.write((33-24-len(crystal_col[3]))*' '+crystal_col[3])
    sys.write((40-33-len(crystal_col[4]))*' '+crystal_col[4])
    sys.write((47-40-len(crystal_col[5]))*' '+crystal_col[5])
    sys.write((54-47-len(crystal_col[6]))*' '+crystal_col[6]+' ')
    sys.write((66-55-len(crystal_col[7]+' '+crystal_col[8]))*' '+crystal_col[7]+' '+crystal_col[8])
    sys.write((70-66-len(crystal_col[9]))*' '+crystal_col[9])
    sys.write('\n')

    for line in sub :
        cols = line.split()
        if cols[0] == "ATOM" :
            sys.write(line)

    n_atom_mol = 3
    n_count = 0

    for line in liq :
        cols = line.split()
        if cols[0] == "ATOM" :
            n_count = ( n_count + 1 ) % n_atom_mol
            if n_count == 1 :
                y = float(cols[6])
                not_carve = y <= sub_max_y
            if not_carve :
                sys.write(line)

    sys.close()
    liq.close()
    sub.close()


"""
    FUNCTION add_atoms_topology
    [...]
"""
def add_atoms_topology(
    system_file = "system.pdb",
    topology_file = "system.top"
    ):

    sys = open(system_file, 'r')

    residues_number = dict()
    types_per_residual = dict()

    for line in sys :
        cols = line.split()
        if cols[0] == "CRYST1" :
            break

    for line in sys :
        cols = line.split()
        residue_name = cols[3]
        if residue_name in residues_number :
            residues_number[residue_name] += 1
        else :
            residues_number[residue_name] = 1
        atom_name = cols[2]
        if residue_name in types_per_residual :
            types_per_residual[residue_name].add(atom_name)
        else :
            types_per_residual[residue_name] = {atom_name}

    sys.close()

    top = open(topology_file, 'a+')

    top.write("[ molecules ]\n")

    for residue_name in residues_number :
        n_mol = residues_number[residue_name]/len(types_per_residual[residue_name])
        top.write(residue_name+' '+str(int(n_mol))+"\n")

    top.close()


"""
    DICTIONARY unit_cell_type
    [...]
"""
unit_cell_type = dict()
unit_cell_type['cube'] = np.array([ [0.0, 0.0, 0.0] ], dtype=np.float64)
unit_cell_type['bcc'] = np.array([ [0.0, 0.0, 0.0], [0.5, 0.5, 0.5] ], dtype=np.float64)
unit_cell_type['fcc'] = np.array([ [0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5] ], dtype=np.float64)


"""
    CLASS monatomic_crystal_substrate
    [...]
"""
class monatomic_crystal_substrate :

    def __init__ ( self, ncx, ncy, ncz, sp, ct='cube' ) :
        self.n_cells_x = ncx
        self.n_cells_y = ncy
        self.n_cells_z = ncz
        self.spacing = sp
        self.crystal_type = ct
        self.local_coordinates = unit_cell_type[ct]
        self.n_atoms_per_cell = (self.local_coordinates.size)//3
        self.n_atoms = self.n_atoms_per_cell*ncx*ncy*ncz
        self.x_coord = np.zeros( self.n_atoms )
        self.y_coord = np.zeros( self.n_atoms )
        self.z_coord = np.zeros( self.n_atoms )
        self.init_coordinates()

    def init_coordinates( self ) :
        print("Initializing "+self.crystal_type+" crystal coordinates ...")
        n = 0
        for i in range(self.n_cells_x) :
            for j in range(self.n_cells_y) :
                for k in range(self.n_cells_z) :
                    for l in range(self.n_atoms_per_cell) :
                        self.x_coord[n] = self.spacing * ( self.local_coordinates[l][0] + i )
                        self.y_coord[n] = self.spacing * ( self.local_coordinates[l][1] + j )
                        self.z_coord[n] = self.spacing * ( self.local_coordinates[l][2] + k )
                        n += 1
        print("Crystal initialized!")

    def carve_function ( self, fun ) :
        x_new = []
        y_new = []
        z_new = []
        for n in range(self.n_atoms) :
            if self.z_coord[n] < fun(self.x_coord[n]) :
                x_new.append(self.x_coord[n])
                y_new.append(self.y_coord[n])
                z_new.append(self.z_coord[n])
        self.x_coord = np.array(x_new)
        self.y_coord = np.array(y_new)
        self.z_coord = np.array(z_new)
        self.n_atoms = self.x_coord.size

    def write_to_file( self, file_name ) :
        fn = open(file_name, 'w')
        n = 1
        fn.write( "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" %
           (self.x_coord.max(),self.y_coord.max(),self.z_coord.max(),90.0,90.0,90.0) )
        for n in range(self.n_atoms):
            a = 'CUS'   # a = 'CUB'
            fn.write( "%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n" %
                ("ATOM",n % 100000,a,"SUB",1,self.x_coord[n],self.y_coord[n],self.z_coord[n]) )
        fn.close()
