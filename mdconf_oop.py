"""
    OOP version of mdconf library
"""

import numpy as np

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
    a = Atom()
    a.mol_idx = int( line[0:5].strip() )
    a.residue = str( line[5:10].strip() )
    a.type = str( line[10:15].strip() )
    a.atom_idx = int( line[15:20].strip() )
    a.x_pos = float( line[20:28].strip() )
    a.y_pos = float( line[28:36].strip() )
    a.z_pos = float( line[36:44].strip() )
    a.x_vel = non_null_float( line[44:52].strip() )
    a.y_vel = non_null_float( line[52:60].strip() )
    a.z_vel = non_null_float( line[60:68].strip() )
    return a

def read_pdb_line( line ) :
    a = Atom()
    a.mol_idx = int( line[22:26].strip() )
    a.residue = str( line[17:20].strip() )
    a.type = str( line[12:16].strip() )
    a.atom_idx = int( line[6:11].strip() )
    a.x_pos = float( line[30:38].strip() )
    a.y_pos = float( line[38:46].strip() )
    a.z_pos = float( line[46:54].strip() )
    # Velocity is not defined in a .pdb file
    a.x_vel = np.nan
    a.y_vel = np.nan
    a.z_vel = np.nan
    return a

def read_box_pdb( file_in, pattern = "CRYST1" ) :
    f = open(file_in, 'r')
    for line in f :
        cols = line.split()
        if cols[0] == pattern:
            bx = float(cols[1])
            by = float(cols[2])
            bz = float(cols[3])
            break
    f.close()
    return bx, by, bz

class MDSystem :

    def __init__ ( MDSystem ) :
        MDSystem.init_file = " "
        MDSystem.tag = "blank_system"
        MDSystem.title = " "
        MDSystem.n_atoms = 0
        MDSystem.atoms = []
        MDSystem.box_xx = 0.0
        MDSystem.box_yy = 0.0
        MDSystem.box_zz = 0.0
        MDSystem.box_xy = 0.0
        MDSystem.box_xz = 0.0
        MDSystem.box_yz = 0.0

    def __init__ ( MDSystem, file_name ) :
        MDSystem.init_file = file_name
        MDSystem.tag = file_name.split('.')[0]
        ext = file_name.split('.')[-1]
        if ext == "gro" :
            print("MDSystem: reading from .gro file")
            MDSystem.read_gro()
        elif ext == "pdb" :
            print("MDSystem: reading from .pdb file (velocity values not included!)")
            MDSystem.read_pdb()
        else :
            raise Exception("Unsupported file type")

    def read_gro ( MDSystem ) :
        # Count number of lines
        n_lines = count_line( MDSystem.init_file )
        # Store content
        f_in = open(MDSystem.init_file, 'r')
        MDSystem.title = f_in.readline()
        MDSystem.n_atoms = int(f_in.readline())
        MDSystem.atoms = []
        for n in range( n_lines-3 ) :
            line = f_in.readline()
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

    def read_pdb ( MDSystem ) :
        print('Warning: initialization via .pdb does not support non-cubic boxes')
        MDSystem.box_xx, MDSystem.box_yy, MDSystem.box_zz = \
            read_box_pdb( MDSystem.init_file )
        MDSystem.box_xx = np.nan
        MDSystem.box_yy = np.nan
        MDSystem.box_zz = np.nan
        MDSystem.atoms = []
        f_in = open(MDSystem.init_file, 'r')
        for line in f_in :
            record_type = line[0:5].strip()
            if record_type=="TITLE" :
                MDSystem.title = line[6:]
            elif record_type=="ATOM" :
                MDSystem.atoms.append( read_pdb_line( line ) )
        MDSystem.n_atoms = len( MDSystem.atoms )
        f_in.close()
