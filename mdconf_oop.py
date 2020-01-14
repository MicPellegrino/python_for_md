"""
    OOP version of mdconf library
"""

class Atom :
    mol_idx = 0
    residue = " "
    type = " "
    atom_idx = 0
    x_pos = 0.0
    y_pos = 0.0
    z_pos = 0.0
    x_vel = 0.0
    y_vel = 0.0
    z_vel = 0.0

def count_line( file_name ) :
    n_lines = 0
    f_in = open(MDSystem.init_file, 'r')
    for line in MDSystem.init_file :
        n_lines += 1
    f_in.close()
    return n_lines

def non_null_float( s ) :
    if s == '' or s == '\n' :
        return = 0.0
    else :
        return float( s )

def read_gro_line( line ) :
    a = Atom()
    a.mol_idx = int(line[0:5])
    a.residue = str(line[5:10])
    a.type = str(line[10:15])
    a.atom_idx = int(line[15:20])
    a.x_pos = float(line[20:28])
    a.y_pos = float(line[28:36])
    a.z_pos = float(line[36:44])
    a.x_vel = non_null_float( line[44:52] )
    a.y_vel = non_null_float( line[52:60] )
    a.z_vel = non_null_float( line[60:68] )
    return a

class MDSystem :

    def __init__ ( MDSystem, file_name ) :
        MDSystem.init_file = file_name
        ext = file_name.split('.')[-1]
        if ext == "gro" :
            MDSystem.read_gro()
        elif ext == "pdb" :
            print("WARNING: pdb files do not contain information regarding atoms velocities")
            MDSystem.read_pdb()
        else :
            raise Exception("Unsupported file type")

    def read_gro ( MDSystem ) :
        # Count number of lines
        n_lines = count_line( file_name )
        # Store content
        f_in = open(MDSystem.init_file, 'r')
        MDSystem.title = f_in.readline()
        MDSystem.n_atoms = int(f_in.readline())
        MDSystem.atoms = []
        for idx in range(3,n_lines) :
            MDSystem.atoms.append( read_gro_line( line ) )
        assert len(MDSystem.atoms) == MDSystem.n_atoms, "ERROR: inconsistent number of atoms in gro file"
        line = f_in.readline()
        MDSystem.box_xx = non_null_float( line[0:10]  )
        MDSystem.box_yy = non_null_float( line[10:20] )
        MDSystem.box_zz = non_null_float( line[20:30] )
        MDSystem.box_xy = non_null_float( line[30:40] )
        MDSystem.box_xz = non_null_float( line[40:50] )
        MDSystem.box_yz = non_null_float( line[50:60] )
        f_in.close()

    def read_pdb () :
        # ...
