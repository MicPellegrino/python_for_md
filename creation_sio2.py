import mdconf as mdc

h = 4.5
h_off = h
w_n = 1.0/h
w_off = 0.0
ni = 100
nj = 20
nk = 1
file_name = 'energy_minimization/quad_substrate_wave.pdb'
bend = True

mdc.quad_substrate_wave(h, h_off, w_n, w_off, bend, ni, nj, nk, file_name)

mdc.add_atoms_topology('energy_minimization/quad_substrate_wave.pdb', 'energy_minimization/quad_substrate_wave.top')
