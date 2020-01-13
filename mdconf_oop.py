"""
    OOP version of mdconf library
"""

class MDSystem :

    def __init__ ( MDSystem, file_name ) :
        MDSystem.init_file = file_name
        ext = file_name.split('.')[-1]
        if ext == "gro" :
            MDSystem.read_gro()
        elif ext == "pdb" :
            MDSystem.read_pdb()
        else :
            # Throw an exception ...
            print("Non so lanciare exceptions con python!")

    def read_gro ( MDSystem ) :
        # Count number of lines
        f_in = open(MDSystem.init_file, 'r')
        n_lines = 0
        for line in MDSystem.init_file :
            n_lines += 1
        f_in.close()
        # Store content
        f_in = open(MDSystem.init_file, 'r')
        idx = 0
        for line in MDSystem.init_file :
            idx += 1
            if idx == 1 :
                MDSystem.title = line
            elif idx == 2 :
                MDSystem.n_atoms = int(line)
            elif idx == n_lines :
                MDSystem.box_xx = float( line[0:10]  )
                MDSystem.box_yy = float( line[10:20] )
                MDSystem.box_zz = float( line[20:30] )
                MDSystem.box_xy = float( line[30:40] )
                MDSystem.box_xz = float( line[40:50] )
                MDSystem.box_yz = float( line[50:60] )
            else :
                # La struttura dati migliore in termini di chiarezza
                # sarebbe una lista di structs
        f_in.close()

    def read_pdb () :
