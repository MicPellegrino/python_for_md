import mdconf_oop as md

root_folder = "/home/michele/python_for_md/PureWater"

input_files = [root_folder+'/DeformEm1/init-conf.gro', 
        root_folder+'/DeformEm2/init-conf.gro', 
        root_folder+'/DeformEm3/init-conf.gro', 
        root_folder+'/DeformEm4/init-conf.gro', 
        root_folder+'/DeformEm5/init-conf.gro']

output_files = [root_folder+'/DeformEm1/shear-conf.gro',
        root_folder+'/DeformEm2/shear-conf.gro', 
        root_folder+'/DeformEm3/shear-conf.gro', 
        root_folder+'/DeformEm4/shear-conf.gro', 
        root_folder+'/DeformEm5/shear-conf.gro']

initial_vel = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

for i in range(len(initial_vel)) :

    print("######################################")
    print("Deformatio rate = "+str(initial_vel[i]))
    conf = md.Configuration()
    conf.input(input_files[i])
    conf.assign_shear_velocity(initial_vel[i])
    z, v = conf.bin_velocity_profile()
    print(v)
    conf.output(output_files[i])
    print("######################################")
