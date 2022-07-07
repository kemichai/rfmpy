import rfmpy.core.migration_sphr as rf_mig
import rfmpy.utils.migration_plots_spher as plot_migration_sphr
import numpy as np
import os
import matplotlib.pyplot as plt
from vtk.util import numpy_support
import vtk


# Define paths
work_dir = os.getcwd()
path = work_dir + "/data/RF/"

with open('obs_amplitudes_matrix.npy', 'rb') as f:
    G3 = np.load(f)



from pyevtk.hl import gridToVTK
import numpy as np
import random as rnd
# # Dimensions
# nx, ny, nz = 6, 6, 2
# lx, ly, lz = 1.0, 1.0, 1.0
# dx, dy, dz = lx/nx, ly/ny, lz/nz
# ncells = nx * ny * nz
# npoints = (nx + 1) * (ny + 1) * (nz + 1)
# # Coordinates
X = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
Y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
Z = np.arange(0, lz + 0.1*dz, dz, dtype='float64')
x = np.zeros((nx + 1, ny + 1, nz + 1))
y = np.zeros((nx + 1, ny + 1, nz + 1))
z = np.zeros((nx + 1, ny + 1, nz + 1))
# We add some random fluctuation to make the grid more interesting
for k in range(nz + 1):
 for j in range(ny + 1):
     for i in range(nx + 1):
         x[i,j,k] = X[i] + (0.5 - rnd.random()) * 0.2 * dx
         y[i,j,k] = Y[j] + (0.5 - rnd.random()) * 0.2 * dy
         z[i,j,k] = Z[k] + (0.5 - rnd.random()) * 0.2 * dz
# # Variables
# pressure = np.random.rand(ncells).reshape( (nx, ny, nz))
# temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1))
gridToVTK("./structured", x, y, z, cellData = {"pressure" : pressure}, pointData = {"temp" : temp})



from pyevtk.hl import imageToVTK

inc = 0.25
zmax = 100
# Determine study area (x -> perpendicular to the profile)
minx = 5.0
maxx = 12.0
pasx = 0.5
miny = 45.0
maxy = 52.0
pasy = 0.5
minz = -2
# maxz needs to be >= zmax
maxz = 100
pasz = 2


x = np.arange(minx, maxx + pasx, pasx)
y = np.arange(miny, maxy + pasy, pasy)
z = np.arange(minz, maxz + pasz, pasz)

# Dimensions
nx, ny, nz = len(x), len(y), len(z)
lx, ly, lz = pasx, pasy, pasz
dx, dy, dz = lx/nx, ly/ny, lz/nz
ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)

# Variables
pressure = G3.reshape( (nx, ny, nz), order = 'C')
temp = G3.reshape( (nx, ny, nz))
gridToVTK("./structured", x, y, z, cellData = {"pressure" : pressure}, pointData = {"temp" : temp})

