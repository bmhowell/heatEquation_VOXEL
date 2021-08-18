"""
Voxel code
Brian Howell
"""

from voxelMESH import *
import matplotlib.pyplot as plt
import scipy.sparse
import scipy.sparse.linalg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import ArtistAnimation
import matplotlib.cm as cm

# constants
thetaW = 300            # |    K    |  temperature at the wall
theta0 = 300            # |    K    |  initial temperature
I0 = 2.0 * 10**8        # |  W/m^2  |  initial laser intensity
a = 25.0                # |   1/m   |  absorption constant
L = 0.05                # |    m    |  sample length
node = 31               # |   ___   |  number of nodes

# material properties
Kp = 10                 # |  W/m-K  |  thermal conductivity
C = 450                 # | J/kg-K  |  heat capacity
rho0 = 3000             # | kg/m^3  |  initial material density
vP = 0.3                # |   ---   |  volume fraction of particles
rPart = L / 10          # |    m    |  radius of the particles

# simulation parameters
tf = 2                  # |    s    |  final simulation time
dt = 1e-4               # |    s    |  initial time step

# compute mesh and A matrix
A, coord = computeCoord(node, L)

print('coord: ', coord)
