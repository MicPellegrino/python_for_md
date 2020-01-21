import numpy as np
import matplotlib.pyplot as plt

# density_array = np.fromfile('densmap.dat', dtype=float)
# N_da = len(density_array)

# size of simulation box
Lx = 47.72970 
Lz = 35.31586
h = 0.05

Nda = 0
for line in open('densmap.dat', 'r'):
    vals = np.array( [ float(i) for i in line.split() ] )
    Nda += len(vals)

print(Nda)

Nx = np.sqrt(Nda*Lx/Lz)
Nz = Nx*Lz/Lx

Nx = int( np.round(Nx) )
Nz = int( np.round(Nz) )

print(Nx*Nz)

idx = 0
density_array = np.zeros(Nx*Nz, dtype=float)
for line in open('densmap.dat', 'r'):
    vals = np.array( [ float(i) for i in line.split() ] )
    density_array[idx:idx+len(vals)] = vals
    idx += len(vals)

# print(Nx*Nz)
# print(len(density_array))

# b = np.sqrt(N_da*Lx/Lz)
# a = (Lz/Lx)*b

# print("a="+str(a))
# print("b="+str(b))
# print("a*b="+str(a*b))
# print("N_da="+str(N_da))

density_array = density_array.reshape((Nx,Nz))
density_array = density_array[int(3*Nx/8):int(5*Nx/8),0:int(Nz/4)]

nx = density_array.shape[0]
nz = density_array.shape[1]

grad_x = np.zeros((nx, nz), dtype=float)
grad_z = np.zeros((nx, nz), dtype=float)

# need only the bulk values
for i in range(1, nx-1) :
    for j in range (1, nz-1) :
        grad_x[i,j] = 0.5*( density_array[i+1,j] - density_array[i-1,j] )/h
        grad_z[i,j] = 0.5*( density_array[i,j+1] - density_array[i,j-1] )/h

grad_norm = np.sqrt( np.power( grad_x, 2 ) + np.power( grad_z, 2 ) )

# Plotting
plt.matshow(np.flip(density_array.transpose(), axis=0))
# plt.contour(density_array.transpose())
plt.show()
plt.matshow(np.flip(grad_norm.transpose(), axis=0))
# plt.contour(grad_norm.transpose())
plt.show()

