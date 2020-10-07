import os
import mdconf as md

main_folder = '/home/michele/python_for_md/Droplet50nmExp'
topology_folder = main_folder+'/Topologies'
topology_prototypes = [ main_folder+'/topology_q1.top',
                        main_folder+'/topology_q2.top',
                        main_folder+'/topology_q3.top',
                        main_folder+'/topology_q4.top',
                        main_folder+'/topology_q5.top'  ]

topology = [ ]

input_height = [    main_folder+'/H1/merge_height1.pdb',
                    main_folder+'/H2/merge_height2.pdb',
                    main_folder+'/H3/merge_height3.pdb',
                    main_folder+'/H4/merge_height4.pdb' ] 
input_lambda = [    main_folder+'/L1/merge_lambda1.pdb',
                    main_folder+'/L2/merge_lambda2.pdb',
                    main_folder+'/L3/merge_lambda3.pdb',
                    main_folder+'/L4/merge_lambda4.pdb' ]

# Create topology files
for q in range(5) :
    for l in range(4) :
        file_name = topology_folder+'/l'+str(l+1)+'q'+str(q+1)+'.top'
        os.system( "cp "+topology_prototypes[q]+" "+file_name )
        topology.append(file_name)
    for h in range(4) :
        file_name = topology_folder+'/h'+str(h+1)+'q'+str(q+1)+'.top'
        os.system( "cp "+topology_prototypes[q]+" "+file_name )
        topology.append(file_name)

# Add atoms to topology files
# Flat

# Rough
n = 0
for q in range(5) :
    for l in range(4) :
        md.add_atoms_topology( input_lambda[l], topology[n] )
        n = n+1
    for h in range(4) :
        md.add_atoms_topology( input_height[h], topology[n] )
        n = n+1

"""
for idx in range(5) :
    file_system = '/home/michele/python_for_md/Droplet20nmExp/Systems/system_wave'+str(int(idx+1))+'.pdb'
    for q in range(5) :
        file_topology = '/home/michele/python_for_md/Droplet20nmExp/Topologies/r'+str(int(idx+1))+'q'+str(int(q+1))+'.top'
        md.add_atoms_topology( file_system, file_topology )
"""
