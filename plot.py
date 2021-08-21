"""
Plot Voxel results
-Brian Howell
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd

# import data
# df = pd.read_csv("voxel3.dat", sep='\t')
df = pd.read_csv("state.dat", sep=' ')
data = np.asarray(df.values)
xPos = data[:, 0]
yPos = data[:, 1]
zPos = data[:, 2]
density = data[:, 3]

sideNode = int(np.round(len(xPos) ** (1 / 3), 0))
x = xPos.reshape((sideNode, sideNode, sideNode))
y = yPos.reshape((sideNode, sideNode, sideNode))
z = zPos.reshape((sideNode, sideNode, sideNode))
d = density.reshape((sideNode, sideNode, sideNode))

# normalize density
dNorm = (d - d.min().min().min())
dNorm = dNorm / dNorm.max().max().max()

# reshape and
xVox = np.ones((x.shape), dtype=bool)
yVox = np.ones((y.shape), dtype=bool)
zVox = np.ones((z.shape), dtype=bool)

voxels = xVox | yVox | zVox
print('voxel: ', voxels.shape)

# set the colors of each object
colors = np.empty(voxels.shape, dtype=object)
ax = plt.figure().add_subplot(projection='3d')
ax.voxels(voxels, facecolors=plt.cm.jet(dNorm), edgecolor='k')

plt.savefig("initParticle.png")
plt.show()

