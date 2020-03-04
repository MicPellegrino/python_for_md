import mdconf as md


# Flat
file_system = '/home/michele/python_for_md/Droplet20nmExp/Systems/system_flat.pdb'
for q in range(5) :
    file_topology = '/home/michele/python_for_md/Droplet20nmExp/Topologies/r0q'+str(int(q+1))+'.top'
    md.add_atoms_topology( file_system, file_topology )

# Rough
for idx in range(5) :
    file_system = '/home/michele/python_for_md/Droplet20nmExp/Systems/system_wave'+str(int(idx+1))+'.pdb'
    for q in range(5) :
        file_topology = '/home/michele/python_for_md/Droplet20nmExp/Topologies/r'+str(int(idx+1))+'q'+str(int(q+1))+'.top'
        md.add_atoms_topology( file_system, file_topology )

