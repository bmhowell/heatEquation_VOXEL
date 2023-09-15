"""
Project 5
ME C201: Modeling and Simulation of Advanced Manufacturing Processes
Spring Semester 2021
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import scipy.sparse.linalg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import ArtistAnimation
import matplotlib.cm as cm
import pandas as pd
import os

# Constants                                    |UNITS    | Description
K = 135  # | W/m-K   | Thermal conductivity
# conductivity
C = 450  # | J/kgK   | heat capacity
rho_0 = 6500  # | kg/m^3  | density
theta_w = 300  # |   K     | temp at the wall
theta_0 = 300  # |   K     | initial temp
I_0 = 2e7  # |  W/m^2  | Initial laser intensity
a = 25  # |  1 / m  | absorption constant
tf = 2  # |   s     | final simulated time
dt = 1e-4  # |   s     | time step size
L = 0.05  # |   mm    | side length of cube

# load data from c++ generation of particles
df = pd.read_csv("density.dat", sep=' ')
data = np.asarray(df.values)
xPos = data[:, 0]
yPos = data[:, 1]
zPos = data[:, 2]
density = data[:, 3]
heatCap = data[:, 4]
thermCond = data[:, 5]

def intensity(zz):
    intense = I_0 * np.exp(-a * (L - zz))
    return intense


def laser(zz):
    energy_output = a * I_0 * np.exp(-a * (L - zz))
    return energy_output


def temp_profile(N):
    # global K, C, rho_0, theta_w, theta_0, I_0, a, tf, dt, L
    h = L / (N - 1)
    print('... Computing A ...')

    # compute coordinates for every node
    coord_cube = np.zeros((N ** 3, 4))  # array containing: Node #, x, y, z
    x = y = z = np.arange(0, L + h, h)
    counter = 0
    for k in range(len(z)):
        for j in range(len(y)):
            for i in range(len(x)):
                coord_cube[counter, 0] = counter
                coord_cube[counter, 1], coord_cube[counter, 2], coord_cube[counter, 3] = x[i], y[j], z[k]
                counter = counter + 1

    top_boundary = []
    top_boundary_h = []
    bottom_boundary = []
    bottom_boundary_h = []
    wall_boundary = []
    counter = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i == 0:
                    bottom_boundary.append(counter)
                    bottom_boundary_h.append(counter + N ** 2)
                if i == N - 1:
                    top_boundary.append(counter)
                    top_boundary_h.append(counter - N ** 2)
                if j == 0 or j == N - 1 or k == 0 or k == N - 1:
                    wall_boundary.append(counter)
                counter += 1

    # Construct A matrix => N^3 x N^3 where N=21
    A = np.zeros((N ** 3, N ** 3))
    for i in range(N ** 3):
        # set components of A (six point stencil) for nodes not on the top or bottom face
        if (N ** 2) < i < ((N ** 3) - (N ** 2)):
            A[i, i] = -6
            A[i, i - 1] = A[i, i + 1] = 1  # capture neighboring nodes in x-direction
            A[i, i - N] = A[i, i + N] = 1  # capture neighboring nodes in y-direction
            A[i, i - N ** 2] = A[i, i + N ** 2] = 1  # capture neighboring nodes in z-direction


    # zero out the rows in A that pertain to nodes on any of the boundaries
    # for i in range(len(boundary_node)):
    #     A[boundary_node[i], :] = 0  # captures Dirichlet BC's. Initializes Neumann BCs
    A[wall_boundary, :] = 0
    A = scipy.sparse.csr_matrix(A)

    print('Finished computing A')

    # Compute the I_abs vector
    I_abs = np.zeros(N ** 3)
    for i in range(len(coord_cube)):
        if 0.02 <= coord_cube[i, 1] <= 0.03 and 0.02 <= coord_cube[i, 2] <= 0.03:
            I_abs[(coord_cube[i, 0]).astype(int)] = laser(coord_cube[i, 3])

    # using forward euler, solve differential equation
    theta = np.ones(N ** 3) * theta_0  # compute initial temperature array:
    theta[wall_boundary] = theta_w
    temp = []

    # r_ims = []
    # fig = plt.figure()

    X, Y = np.meshgrid(x, y)
    t = np.arange(0, tf + dt, dt)  # discretize time
    for i in range(len(t)):
        print('average temp: ', np.average(theta))

        # for some reason the following line does not work
        # theta = theta + dt / rho_0 / C * (K / h ** 2 * A @ theta + I_abs)

        Atheta = A @ theta  # but if you precompute this term it will work
        theta = theta + dt / density / heatCap * (thermCond / h ** 2 * Atheta + I_abs)

        # enforce Neumman BCs
        theta[bottom_boundary] = theta[bottom_boundary_h]
        theta[top_boundary] = theta[top_boundary_h]

        # for mp4
        # theta_plot = []
        # for i in range(N ** 3 - N ** 2, N ** 3):
        #     theta_plot.append(theta[i])
        # theta_plot = np.reshape(theta_plot, (N, N))
        # cp = plt.contourf(X, Y, theta_plot, cmap='hot')
        # r_ims.append(cp.collections)

        if i % 100 == 0:
            temp.append(theta)

    # create mp4
    # cbar = plt.colorbar(cp)
    # cbar.set_label(r'$\theta$ (K)')
    # plt.axes().set_aspect('equal')
    # ani = ArtistAnimation(fig, r_ims, interval=1, repeat=True)
    # ani.save('temperature_time.mp4', writer='ffmpeg')

    # plot average of top layer for every time step
    average_top = []
    average_bottom = []
    for i in range(len(temp)):
        average_top.append(np.average(temp[i][top_boundary]))
        average_bottom.append(np.average(temp[i][bottom_boundary]))

    # print('average_top = ', average_top)
    # print('average_bottom = ', average_bottom)

    time_steps = np.arange(0, len(temp), 1) / 1000
    plt.figure()
    plt.plot(time_steps, average_bottom, label='Bottom face')
    plt.plot(time_steps, average_top, label='Top face')
    plt.legend()
    plt.title('Average Temperatures: Top & Bottom Faces')
    plt.ylabel(r'$\theta$ (K)')
    plt.xlabel(r'Time (s)')
    plt.savefig('average.png')

    # plot contour of top layer of cube
    top_cube_temp = []
    for i in range(N ** 3 - N ** 2, N ** 3):
        top_cube_temp.append(temp[-1][i])

    top_temp = np.reshape(top_cube_temp, (N, N))
    plt.figure()
    cp = plt.contourf(X, Y, top_temp, 15, cmap=cm.hot, vmin=np.amin(top_temp), vmax=np.amax(top_temp))
    m = plt.cm.ScalarMappable(cmap=cm.hot)
    m.set_array(top_temp)
    m.set_clim(np.amin(top_temp), np.amax(top_temp))
    plt.ylabel('y(m)')
    plt.xlabel('x(m)')
    cbar = plt.colorbar(m, boundaries=np.linspace(np.amin(top_temp), np.amax(top_temp), 100))
    cbar.set_label("Temperature (K)")

    plt.title(r'$\theta$ (z = 0.05) distribution on XY plane')
    plt.savefig('laser_output_top.png')

    # plot contour of top layer of cube
    bottom_cube_temp = []
    for i in range(0, N ** 2):
        bottom_cube_temp.append(temp[-1][i])

    bottom_temp = np.reshape(bottom_cube_temp, (N, N))


    plt.figure()
    cp = plt.contourf(X, Y, bottom_temp, 15, cmap=cm.hot, vmin=np.amin(top_temp), vmax=np.amax(top_temp))
    m = plt.cm.ScalarMappable(cmap=cm.hot)
    m.set_array(top_temp)
    m.set_clim(np.amin(top_temp), np.amax(top_temp))
    plt.ylabel('y(m)')
    plt.xlabel('x(m)')
    cbar2 = plt.colorbar(m, boundaries=np.linspace(np.amin(top_temp), np.amax(top_temp), 100))
    cbar2.set_label("Temperature (K)")

    plt.title(r'$\theta$ (z = 0.00) distribution on XY plane')
    plt.savefig('laser_output_bottom.png')

    return temp, coord_cube


tempTotal, cubeCoord = temp_profile(31)

# #save voxels - full grid density comparison
# if os.path.exists("stateTemp.dat"):
#     os.remove("stateTemp.dat")
#     # Print the statement once
#     # the file is deleted
#     print("Previous file deleted.")
save_path1 = "/Users/bhopro/Desktop/Berkeley/MSOL/Projects/Voxel/stateTemp.dat"

f = open(save_path1, "w+")
f.write('X Y Z T \n')
# for i in range(len(tempTotal)):
#     f.write(', T{}'.format(i+1))
# f.write('\n')

for i in range(len(cubeCoord[:, 0])):
    f.write("{} {} {} {} \n".format(cubeCoord[i, 1], cubeCoord[i, 2], cubeCoord[i, 3], tempTotal[-1][i]))




f.close()




