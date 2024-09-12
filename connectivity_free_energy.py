import mdconf_oop as md
import numpy as np

# Volume (should stay constant since we simulate a NVT ensemble)
V0 = 18.59721*4.67654*7.50000

R0 = 25.0
r_frac = 0.15
n = 5

l = R0*r_frac
k = 2.0*np.pi/l

Ly = 4.67654
a = 0.7
Lx = n*l
Lz = V0/(Lx*Ly)
conf = md.Configuration(Lx, Ly, Lz)
conf.silica_monolayer_rough(0.5*Lz, a/k, k, 0.0, cut_box=True)
conf.print_info()
conf.output('connectR15.gro')

nj = int(Ly/(0.450*md.alpha_1))
ni = int((conf.n_atoms/3.0)/nj)
print(nj*ni)

Lx = ni*0.450
Lz = V0/(Lx*Ly)
conf = md.Configuration(Lx, Ly, Lz)
conf.silica_monolayer(0.5*Lz)
conf.print_info()
conf.output('connectFL.gro')
