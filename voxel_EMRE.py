'''Voxel Code for Thermodynamics Simulation.
Code written by: Emre Mengi.
August 2021.'''

################################## IMPORTING PACKAGES ####################################################

import numpy as np
import os
import math
import csv
from scipy import sparse

################################## FUNCTIONS ##############################################################

######################################## Mesh Initialization and Material Creation ########################

def myCoord(Nmesh,L): #create coordinate array for mesh

    x = np.linspace(0,L,Nmesh) #nodal coordinate in 1-D
    Xplot, Yplot = np.meshgrid(x,x) #meshgrid for plotting
    coord = np.zeros([Nmesh**3,3]); #pre-allocate coordinate matrix
    for j in range(Nmesh):
        for i in range(Nmesh**2):
            coord[(Nmesh**2)*j+i,0] = Yplot[np.unravel_index(i, Yplot.shape, 'F')]
            coord[(Nmesh**2)*j+i,1] = Xplot[np.unravel_index(i, Xplot.shape, 'F')]
    adv = 0

    for i in range(Nmesh**3): #loop to populate z-coordinates
        if i <= Nmesh**2-1:
            coord[i,2] = x[adv]
        elif (i % Nmesh**2) == 0:
            adv = adv + 1
            coord[i:i+Nmesh**2,2] = x[adv]

    # for i in range(Nmesh**2,Nmesh**3):
    #     coord[i,2] = coord[i,2] + x[1]-x[0];

    return (coord)

def myParticle(r_par,v_p,h,Nmesh,coord,rho_0m,rho_0p,Cm,Cp,Km,Kp,int_t): #create the material

    #r_node_par = round(r_par / h) #rounded nodal radius [nodal units]
    N_node_par = round(Nmesh**3 * v_p) #number of particle nodes needed for specified volume fraction
    #N_node_par = N_node_par.astype(int)
    #V_par = 4/3 * np.pi * r_node_par**3 #volume per particle
    #N_par = round(N_node_par / V_par) #number of particles

    #density,heat capacity,conductivity map
    dense = np.ones([Nmesh**3,1]) * rho_0m;
    C = np.ones([Nmesh**3,1]) * Cm
    K = np.ones([Nmesh**3,1]) * Km

    if v_p != 0:
        # create 3D particles
        par_ind_list = []
        intermed_ind_list = []
        par_ind = 0
        i = 0
        while np.size(par_ind) < N_node_par and i < 200:
            par_ind_temp = np.ones([Nmesh**3,1])*(-1)
            intermed_ind_temp = np.ones([Nmesh**3,1])*(-1)
            k = 0
            l = 0
            print('Generating particle ',i+1)
            ind_center = math.ceil(np.random.rand()*Nmesh**3)
            for j in range(Nmesh**3):
                if np.sqrt((coord[j][0] - coord[ind_center][0])**2 + (coord[j][1] - coord[ind_center][1])**2 + (coord[j][2] - coord[ind_center][2])**2) <= r_par:
                    par_ind_temp[k] = j
                    k = k + 1
                elif int_t != 0 and np.sqrt((coord[j][0] - coord[ind_center][0])**2 + (coord[j][1] - coord[ind_center][1])**2 + (coord[j][2] - coord[ind_center][2])**2) <= (r_par + int_t*h):
                    intermed_ind_temp[l] = j
                    l = l + 1



            par_ind_list.append(par_ind_temp[:k])
            par_ind = np.vstack(par_ind_list)
            par_ind = np.unique(par_ind)
            intermed_ind_list.append(intermed_ind_temp[:l])
            i = i + 1


        par_ind = np.vstack(par_ind_list)
        par_ind = np.unique(par_ind)
        par_ind = par_ind.astype(int)

        intermed_ind = np.vstack(intermed_ind_list)
        intermed_ind = np.unique(intermed_ind)
        intermed_ind = intermed_ind.astype(int)

        #density,heat capacity,conductivity map
        dense[par_ind] = rho_0p
        C[par_ind] = Cp
        K[par_ind] = Kp

        #adjust matl properties for the particle indices
        if int_t != 0:
            intermed_ind = np.array(list(set(intermed_ind) - set(par_ind))) #remove interface indices that coincide with other particle core indices
            dense[intermed_ind] = (rho_0p + rho_0m) / 2
            K[intermed_ind] = (Kp + Km) / 2
            C[intermed_ind] = (Cp + Cm) / 2




    return (dense,K,C,par_ind)

##################################### Solving Heat Equation #############################################

def myLaserAbsCreate(I_0,a,z_top,x,y,z): #calculate laser absorption values

    if x >= 0.02 and x <= 0.03 and y >= 0.02 and y <= 0.03:
        I_abs = a * I_0 * math.exp(-a*(z_top - z))
    else:
        I_abs = 0

    return (I_abs)

def myLaserAbsAssign(I_0,a,z_top,coord): #assign laser absorption values

    # Assign I_abs values
    I_abs = np.empty([Nmesh**3,1]) #pre-allocate space
    for i in range(Nmesh**3):
        I_abs[i] = myLaserAbsCreate(I_0,a,z_top,coord[i,0],coord[i,1],coord[i,2]) #laser absorption

    return (I_abs)

def myHeatBoundaries(coord): #set BCs for heat eqn

    # Find BC locations
    DBCind = np.hstack([np.where(coord[:,0]==0),np.where(coord[:,0]==0.05),np.where(coord[:,1]==0)\
                      ,np.where(coord[:,1]==0.05)]) #find locations of Dirichlet BC locations
    DBCind = np.unique(DBCind) #remove repeating variables

    NBCind = np.hstack([np.where(coord[:,2]==0),np.where(coord[:,2]==0.05)]); #find locations of Neumann BC locations
    NBCind = np.unique(NBCind); #remove repeating variables

    return (DBCind,NBCind)

def myHeatStencil(Nmesh,DBCind): #set A matrix for heat eqn

    # Create A matrix - DOUBLE CHECK VALUES
    A = np.zeros([Nmesh**3,Nmesh**3]);
    for i in range(Nmesh**3):
        if i >= Nmesh**2:
            A[i,i-Nmesh**2] = 1

        if i >= Nmesh:
            A[i,i-Nmesh] = 1

        if i >= 1:
            A[i,i-1] = 1

        A[i,i] = -6

        if i < Nmesh**3 - Nmesh**2:
            A[i,i+Nmesh**2] = 1

        if i < Nmesh**3 - Nmesh:
            A[i,i+Nmesh] = 1

        if i < Nmesh**3 - 1:
            A[i,i+1] = 1

    A[DBCind,:] = 0 #zero out rows assoc. with side walls
    A = sparse.csr_matrix(A)

    return (A)

def myHeatTimeStep(coord,theta,del_t,dense,C,h,K,A,I_abs,NBCind,Nmesh): #forward euler time stepping - single step


    RHS = del_t / (dense * C) * ((K / h**2) * (A @ theta) + I_abs)
    theta_new = theta + RHS
    for i in range(np.size(NBCind)):
        if NBCind[i]>= Nmesh**2:
            theta_new[NBCind[i]] = theta_new[NBCind[i]-Nmesh**2]
        else:
            theta_new[NBCind[i]] = theta_new[NBCind[i]+Nmesh**2]

        # ind1 = np.where(coord[:,0]==0);
        # ind2 = np.where(coord[:,0]==0.05);
        # theta_new[ind1] = 3000;
        # theta_new[ind2] = 3000;

    return (theta_new)

def myHeatSolverExplicit(theta_0,Nmesh,A,t_f,del_t): #time-stepping for heat eqn - forward euler(explicit)

    theta = theta_0 * np.ones([Nmesh**3,1]) #initial temp array


    #save voxels - full grid density comparison
    if os.path.exists("voxel5.dat"):
        os.remove("voxel5.dat")
        # Print the statement once
        # the file is deleted
        print("Previous file deleted.")

    ## Time-stepping
    t = np.arange(0,t_f,del_t) #time stepping
    N_t = np.size(t) #number of steps
    frame_mod = math.floor(N_t/60) #number of frames to be saved to .dat (60 fps)
    for i in range(N_t):
        print('Time step ',i,' of ',N_t)
        theta = myHeatTimeStep(coord,theta,del_t,dense,C,h,K,A,I_abs,NBCind,Nmesh) #forward-euler time stepping
        if np.mod(i,frame_mod)==0:
            with open('voxel5.dat', 'a') as p_file:
                writer = csv.writer(p_file, delimiter='\t')  # make write variable
                if i == 0:
                    writer.writerow(['X-coord', 'Y-coord', 'Z-coord', 'Temperature','Density'])  # headings
                writer.writerow(['ZONE'])  # self explanatory
                for p in range(Nmesh**3):  # this for-loop writes each particle's (row) data
                    writer.writerow([coord[p,0],coord[p,1],coord[p,2],theta[p,0],dense[p,0]])


    return (theta)

def myHeatTimeStepCorrection(theta,theta_0,L,K): #calculate required time step size correction - CURRENTLY DOESNT WORK

    del_t_correct = 0.1 * 1 / (max(K) * L**2 * (max(theta) - theta_0) / L)

    return(del_t_correct)

############################################################################################################

## Load variables

theta_w = 300 #wall temperature [K]
theta_0 = 300 #initial temperature [K]
I_0 = 2.0 * 10**8 #initial laser intensity [W/m2]
#I_0 = 0 #initial laser intensity [W/m2]
a = 25.0 #absorption constant [1/m]
L = 0.05 #material length [m]
z_top = L #top of z-direction [m]
Nmesh = 31 #mesh size in one dimension
h = L / (Nmesh - 1) #mesh resolution [m]

#particle material properties
Kp = 10 #thermal conductivity (isotropic) [W/m-K]
Cp = 5 #heat capacity [J/kg-K]
rho_0p = 6500 #density [kg/m3]
v_p = 0.3 #approximate volume fraction of particles
r_par = L / 10 #radius of the particles [m]



#N_par = 10 #number of particles

#matrix material properties
Km = 135 #thermal conductivity (isotropic) [W/m-K]
Cm = 450 #heat capacity [J/kg-K]
rho_0m = 3000 #density [kg/m3]

int_t = 1 #interface thickness [nodal units]

t_f = 2 #final time for simulation
del_t = 1e-4 #initial time step size

# Create material and BCs
coord = myCoord(Nmesh,L) #coordinate array
aljklj = klaj
dense,K,C,par_ind = myParticle(r_par,v_p,h,Nmesh,coord,rho_0m,rho_0p,Cm,Cp,Km,Kp,int_t) #initialize particles/matrix



#save voxels - full grid density comparison
if os.path.exists("voxel3.dat"):
    os.remove("voxel3.dat")
    # Print the statement once
    # the file is deleted
    print("Previous file deleted.")
with open('voxel3.dat', 'a') as p_file:
    writer = csv.writer(p_file, delimiter='\t')  # make write variable
    #if j == 0:
    writer.writerow(['X-coord', 'Y-coord', 'Z-coord', 'Density'])  # headings
    writer.writerow(['ZONE'])  # self explanatory
    for p in range(Nmesh**3):  # this for-loop writes each particle's (row) data
      writer.writerow([coord[p,0],coord[p,1],coord[p,2],dense[p][0]])

# # Solve Heat Equation
# (DBCind,NBCind) = myHeatBoundaries(coord) #set BCs
# I_abs = myLaserAbsAssign(I_0,a,z_top,coord) #call laser calc function
# A = myHeatStencil(Nmesh,DBCind) #create A matrix / stencil
# theta = myHeatSolverExplicit(theta_0,Nmesh,A,t_f,del_t) #solve heat equation

# del_t_correct = myHeatTimeStepCorrection(theta,theta_0,L,K)

# # Solve Heat Equation again with Correct Time Step
# I_abs = myLaserAbsAssign(I_0,a,z_top,coord) #call laser calc function
# A = myHeatStencil(Nmesh,DBCind) #create A matrix / stencil
# theta = myHeatSolverExplicit(theta_0,Nmesh,A,t_f,del_t_correct) #solve heat equation
