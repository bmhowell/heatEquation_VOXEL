import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import time

def computeCoord(mSize, L):
    """
    Computes a mesh for the cube of material
    :param mSize: number of nodes along each x-y-z axis
    :return:
    """
    h = L / (mSize - 1)                             # |   ___   |  mesh resolution

    print('---- computing A matrix ----')
    startTime = time.perf_counter()
    coordCube = np.zeros((mSize**3, 4))             # array containing: Node #, x, y, z
    x = y = z = np.arange(0, L+h, h)

    # populate array coordCube with x-y-z coordinates
    counter = 0
    for i in range(len(z)):
        for j in range(len(y)):
            for k in range(len(x)):
                coordCube[counter, 0] = counter
                coordCube[counter, 1] = x[k]
                coordCube[counter, 2] = y[j]
                coordCube[counter, 3] = z[i]
                counter += 1

    # find boundary nodes
    topBoundary = []
    topBoundary_h = []
    bottomBoundary = []
    bottomBoundary_h = []
    wallBoundary = []
    counter = 0
    for i in range(mSize):
        for j in range(mSize):
            for k in range(mSize):
                if i == 0:
                    bottomBoundary.append(counter)
                    bottomBoundary_h.append(counter + mSize**2)
                if i == mSize - 1:
                    topBoundary.append(counter)
                    topBoundary_h.append(counter + mSize**2)
                if j == 0 or j == mSize - 1 or k == 0 or k == mSize - 1:
                    wallBoundary.append(counter)
                counter += 1

    # construct A matrix
    A = np.zeros((mSize**3, mSize**3))
    for i in range(mSize**3):
        # set components of A (six point stencil for nodes not on the top or bottom face
        if (mSize**2) < i < ((mSize**3) - (mSize**2)):
            # avoids top face
            A[i, i] = -6
            A[i, i - 1] = A[i, i + 1] = 1       # capture neighboring nodes in the x-direction
            A[i, i - mSize] = 1
            A[i, i + mSize] = 1                 # capture neighboring nodes in the y-direction
            A[i, i - mSize**2] = 1              # capture neighboring nodes in teh z-direction
            A[i, i + mSize**2] = 1
    # zero out the rows in A that pertain to nodes on any of the boundaries
    A[wallBoundary, :] = 0

    A = scipy.sparse.csr_matrix(A)
    print('---- finished computing A ----')
    print('---- compute time: {} ----'.format(time.perf_counter() - startTime))

    return A, coordCube[:, 1:]



# test function
if __name__ == "__main__":
    from voxelMain import *

    AA, test = computeCoord(21, L)
    print(AA.shape)