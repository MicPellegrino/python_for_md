import mdconf as md

work_dir = '/home/michele/python_for_md/Droplet20nmExp/'

# Flat
f_sys = work_dir+'Equilibration/eq_flat.gro'
f_sh = work_dir+'Systems/sh_flat.gro'
f_vl = work_dir+'Systems/vl_flat.gro'
md.shift_droplet_gro(-6.0, 'z', f_sys, f_sh)
md.add_velocity_gro('z', -0.01, f_sh, f_vl)

# Rough
for idx in range(5) :
    f_sys = work_dir+'Equilibration/eq_wave'+str(idx+1)+'.gro'
    f_sh = work_dir+'Systems/sh_wave'+str(idx+1)+'.gro'
    f_vl = work_dir+'Systems/vl_wave'+str(idx+1)+'.gro'
    md.shift_droplet_gro(-6.0, 'z', f_sys, f_sh)
    md.add_velocity_gro('z', -0.01, f_sh, f_vl)

