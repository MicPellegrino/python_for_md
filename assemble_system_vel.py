import mdconf as md

vel = -0.005

main_folder = '/home/michele/python_for_md/Droplet50nmExp'

input_files = [     main_folder+'/H1/equil_h1.gro',
                    main_folder+'/H2/equil_h2.gro',
                    main_folder+'/H3/equil_h3.gro',
                    main_folder+'/H4/equil_h4.gro'  ]

output_files = [    main_folder+'/H1/veloc_h1.gro',
                    main_folder+'/H2/veloc_h2.gro',
                    main_folder+'/H3/veloc_h3.gro',
                    main_folder+'/H4/veloc_h4.gro'  ]

for h in range(4) :
    md.add_velocity_gro( 'z', vel, input_files[h], output_files[h] )

