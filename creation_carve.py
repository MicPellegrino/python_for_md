## PYTHON FOR MD ##
"""
    Script for the creation of a simple LJ substrate carved according to a 
    fractal structure function
"""

import mdconf as mdc
from math import sqrt

import roughfrac as rf
import numpy as np
from math import log

# Instantiation of 'rough' function
G = 0.5
D = 1.52
gamma = 1.5
x_0 = 100
h_0 = 30
scale_factor = 100
L_micro = G
L_macro = scale_factor*G
n_0 = round( -log(L_macro)/log(gamma) )
n_inf = round( -log(L_micro)/log(gamma) )
test_fun = rf.weierstrass_mandelbrot(G, D, gamma, n_0, n_inf, x_0, h_0)

# Instantiation of crystal structure
ncx = 200
ncy = 10
ncz = 20
sp = 2.7*sqrt(2)
crystal_type = 'fcc'
lj_crystal_struct = mdc.monatomic_crystal_substrate(ncx, ncy, ncz, sp, crystal_type)

# f = lambda x : 0.2*x
# f = lambda x : 40 + 35*(round(x/20.0)%2)
f = lambda x : test_fun.evaluate_1D(x)
lj_crystal_struct.carve_function(f)

file_name = 'lj_sub.pdb'
lj_crystal_struct.write_to_file(file_name)
