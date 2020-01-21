import numpy as np
import matplotlib.pyplot as plt

# density_array = np.fromfile('densmap.dat', dtype=float)
# N_da = len(density_array)

# size of simulation box
Lx = 47.72970 
Lz = 35.31586
h = 0.01
Nx = int( np.round(Lx/h) )
Nz = int( np.round(Lz/h) )

density_array = np.zeros(Nx*Nz, dtype=float)
idx = 0
for line in open('densmap.dat', 'r'):
    vals = np.array( [ float(i) for i in line.split() ] )
    density_array[idx:idx+len(vals)] = vals
    idx += len(vals)

print(Nx*Nz)
print(len(density_array))

# b = np.sqrt(N_da*Lx/Lz)
# a = (Lz/Lx)*b

# print("a="+str(a))
# print("b="+str(b))
# print("a*b="+str(a*b))
# print("N_da="+str(N_da))

density_array = density_array.reshape((Nx,Nz))
density_array = density_array[int(3*Nx/8):int(5*Nx/8),0:int(Nz/4)]
plt.matshow(np.flip(density_array.transpose(), axis=0))
plt.show()
