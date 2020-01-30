import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate
from scipy.interpolate import CubicSpline
import math

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

# SOLVED!
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

# SOLVED!
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

# Set-up values
half_den = 50.0
epsilon = 3.5
rectification = 1000
smoothing_tune = 1e-5

in_range = lambda x: (x<(half_den+epsilon))*(x>(half_den-epsilon))
f = np.vectorize(in_range)
cont_points = f(density_array)
x_data = []
z_data = []
f_data = []
for i in range(0,nx):
    for j in range (0,nz):
        if cont_points[i,j] == 1 :
            x_data.append(h*i)
            z_data.append(h*j)
            f_data.append(density_array[i,j])

m = len(x_data)
x_data_ord = [x_data[0]]
x_data.pop(0)
z_data_ord = [z_data[0]]
z_data.pop(0)
f_data_ord = [f_data[0]]
f_data[0]
dist = np.sqrt( np.power(np.array(x_data)-x_data_ord[0], 2)
    + np.power(np.array(z_data)-z_data_ord[0], 2) )
l = np.argmin(dist)
x_data_ord.append(x_data[l])
x_data.pop(l)
z_data_ord.append(z_data[l])
z_data.pop(l)
f_data_ord.append(f_data[l])
f_data.pop(l)
for k in range(2,m) :
    dir_x_prev = x_data_ord[k-1]-x_data_ord[k-2]
    dir_z_prev = z_data_ord[k-1]-z_data_ord[k-2]
    dist_x = np.array(x_data)-x_data_ord[k-1]
    dist_z = np.array(z_data)-z_data_ord[k-1]
    p = ( ( dir_x_prev*dist_x + dir_z_prev*dist_z ) < 0 ).astype(int)
    dist = np.sqrt( np.power(dist_x, 2) + np.power(dist_z, 2) ) + rectification*p
    l = np.argmin(dist)
    x_data_ord.append(x_data[l])
    x_data.pop(l)
    z_data_ord.append(z_data[l])
    z_data.pop(l)
    f_data_ord.append(f_data[l])
    f_data.pop(l)
x_data_ord.append(x_data_ord[0])
z_data_ord.append(z_data_ord[0])
f_data_ord.append(f_data_ord[0])
x_data_ord = np.array(x_data_ord)
z_data_ord = np.array(z_data_ord)
f_data_ord = np.array(f_data_ord)

weight = -np.absolute(f_data_ord-half_den)/epsilon + 1.0
weight = weight / np.sum(weight)

# This approach may create aliasing at the initial (that is also the final) point
tck, u = interpolate.splprep([x_data_ord, z_data_ord], w=weight, s=smoothing_tune, task=0)
unew = np.arange(0, 1.0, 0.005)
out = interpolate.splev(unew, tck)
x_spline = out[0]
z_spline = out[1]
# np.append(x_spline, x_spline[0])
# np.append(z_spline, z_spline[0])

##############
## Plotting ##
##############

# plt.matshow(np.flip(cont_points.transpose(), axis=0), cmap=cm.bone)
# plt.show()

plt.pcolor(X, Z, density_array, cmap=cm.bone)
plt.colorbar()
plt.plot(x_data_ord, z_data_ord, 'rx', markersize=5)
plt.plot(x_spline, z_spline, 'r-')
# plt.contour(np.flip(density_array.transpose(), axis=0))
plt.axis('scaled')
plt.show()

# plt.hist(density_array.reshape((nx*nz,1)), bins=int(np.sqrt(nx*nz)))
# plt.show()

# plt.matshow(np.flip(grad_norm.transpose(), axis=0), cmap=cm.bone)
# plt.contour(grad_norm.transpose())
# plt.colorbar()
# plt.show()
