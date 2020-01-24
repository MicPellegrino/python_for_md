import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
from scipy.interpolate import CubicSpline

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

Nx = int( np.round(Lx/h) )
Nz = int( np.round(Lz/h) )

# Nx = np.sqrt(Nda*Lx/Lz)
# Nz = Nx*Lz/Lx
# Nx = int( np.round(Nx) )
# Nz = int( np.round(Nz) )

idx = 0
density_array = np.zeros( (Nx+1) * (Nz+1), dtype=float)
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

density_array = density_array.reshape((Nx+1,Nz+1))
density_array = density_array[1:-1,1:-1]
density_array = density_array[int(3*Nx/8):int(5*Nx/8),0:int(Nz/4)]

nx = density_array.shape[0]
nz = density_array.shape[1]

# Meshgrid
x = h*np.arange(0.0,nx,1.0, dtype=float)
z = h*np.arange(0.0,nz,1.0, dtype=float)
X, Z = np.meshgrid(x, z, sparse=False, indexing='ij')

grad_x = np.zeros((nx, nz), dtype=float)
grad_z = np.zeros((nx, nz), dtype=float)

# need only the bulk values
for i in range(1, nx-1) :
    for j in range (1, nz-1) :
        grad_x[i,j] = 0.5*( density_array[i+1,j] - density_array[i-1,j] )/h
        grad_z[i,j] = 0.5*( density_array[i,j+1] - density_array[i,j-1] )/h

grad_norm = np.sqrt( np.power( grad_x, 2 ) + np.power( grad_z, 2 ) )

# Test for spline interpolation
half_den = 50.0
epsilon = 5.0
in_range = lambda x: (x<(half_den+epsilon))*(x>(half_den-epsilon))
f = np.vectorize(in_range)
cont_points = f(density_array)
x_data = []
z_data = []
for i in range(0,nx):
    for j in range (0,nz):
        if cont_points[i,j] == 1 :
            x_data.append(i)
            z_data.append(j)
x_data = h*np.array(x_data)
z_data = h*np.array(z_data)
# tck, u = interpolate.splprep([x_data, z_data], s=0)
# unew = np.arange(0, 1.01, 0.01)
# out = interpolate.splev(unew, tck)

##############
## Plotting ##
##############

# plt.matshow(np.flip(cont_points.transpose(), axis=0), cmap=cm.bone)
# plt.show()

plt.pcolor(X, Z, density_array, cmap=cm.bone)
plt.colorbar()
plt.plot(x_data, z_data, 'rx', markersize=3.5)
# plt.contour(np.flip(density_array.transpose(), axis=0))
plt.show()

# plt.hist(density_array.reshape((nx*nz,1)), bins=int(np.sqrt(nx*nz)))
# plt.show()

# plt.matshow(np.flip(grad_norm.transpose(), axis=0), cmap=cm.bone)
# plt.contour(grad_norm.transpose())
# plt.colorbar()
# plt.show()
