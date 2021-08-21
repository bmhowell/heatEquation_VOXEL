//
// Voxel code CPP
// Created by Brian Howell on 8/17/21.
//


#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <random>
#include <vector>
#include <algorithm>

// particle parameters
const float THETAW = 300.;                              // |    K    |  temperature at the wall
const float THETA0 = 300.;                              // |    K    |  initial temperature
const double I0 = 2 * pow(10, 8);         // |  W/m^2  |  initial laser intensity
const float ABSORB = 25.0;                              // |   1/m   |  absorption constant
const float LENGTH = 0.05;                              // |    m    |  sample length
const int NODE = 7;                                    // |   ___   |  number of nodes
const double INTTHICK = 1.0;                            // |   ---   |  interfacial thickness parameter

// material parameters
const int KP = 10;                                      // |  W/m-K  |  thermal conductivity of particle
const int KM = 135;                                     // |  W/m-K  |  thermal conductivity of material
const int CP = 5;                                       // | J/kg-K  |  heat capacity of particle
const int CM = 450;                                     // | J/kg-K  |  heat capacity
const int RHO0P = 6500;                                 // | kg/m^3  |  initial particle density
const int RHO0M = 3000;                                 // | kg/m^3  |  initial material density
const double VP = 0.5;                                  // |   ---   |  volume fraction of particles
const double RPART = LENGTH / 10.;                      // |    m    |  radius of the particles

// simulation parameters
const double TF = 2.;                                   // |    s    |  final simulation time
const double DT = pow(10, -4);            // |    s    |  initial time step
int SIZETIME = (int) TF / DT;                           // |   ---   |  size of time array
double H = LENGTH / (NODE - 1);                         // physical discretization length
const int SIZEA3 = (int) pow(NODE, 3);           // try static_cast
const int SIZEA2 = (int) pow(NODE, 2);

// data outputs
std::ofstream printState;

// function declarations
void computeCoord(std::vector< std::vector<double> >&);
// computeCoord - compute the coords in the x-y-z directions
// @param vector< vector<double> > - the 2D vector to be filled
// @param double N -  Number of nodes
// @param double L -  Sample length

void computeBoundary(std::vector< std::vector<int> >&);
// boundaryNodes - computes which nodes are on the sides, top and bottom of the cube
// @params vector< vector<double> > - input 2D vector to be filled
// @params double - number of nodes
// @return vector< vector<double> > - vector of vectors containing each side

void computeAsparse(std::vector< std::vector<int> >&,
                    std::vector< std::vector<int> >&);
// computeAsparse - computes the A matrix, in the sparse variety (see comment below for structure)
// @param - modifies 2D vector into a [n x 4] vector

void computeParticles(std::vector< std::vector<double> >&,
                      std::vector<int>&,
                      std::vector<int>&,
                      double[SIZEA3], double[SIZEA3], double[SIZEA3]);
// computeParticles - compute random particles within the medium of material
// @param - modifies 2D vector and filles with particles

void solutionScheme(std::vector< std::vector<int> >&,
                    std::vector< std::vector<int> >&,
                    std::vector< std::vector<double> >&,
                    std::vector<double>&,
                    double [SIZEA3], double [SIZEA3],
                    double [SIZEA3], double [SIZEA3]);
// solutionScheme - iteratitive solver in time and space
// @param - ASparse: FDM mesh matrix
// @param - bNodes: nodes of the boundaries
// @param - temperature: solution vector for every time step
// @param - theta: solution vector for a single time step
// @param - density: array for density of each node
// @param - heatCap: array for heat capacity of each node
// @param - thermCond: array for thermal conductivity
// @param - laserValues: energy input from laser at each node

void write2file(std::vector< std::vector<double> >&,
                double density[SIZEA3]);
// write2file - write cube coordinates and density to file for plotting
// @param - cube coordinates [nParticles x 4] -> [node, x, y, z]
// @param - density [nParticles x 1]

void uniqueVec(std::vector<int> &vec);
// uniqueVec - modifies input vector to only contain unique indices

void print2dVecInt(std::vector< std::vector<int> >&);
// print2dVec - prints 2D vector

void print2dVecDouble(std::vector< std::vector<double> >&);
// print2dVec - prints 2D vector

int main(){
    // initialize time
    double time[SIZETIME];
    double counter = 0.0;
    for (int i = 0; i < SIZETIME; i++) {time[i] = counter; counter += DT;}

    // initialize cubeCoord, bNode, and A matrix vectors
    std::vector< std::vector<double> > cubeCoord;
    std::vector< std::vector<int> > bNodes;
    std::vector< std::vector<int> > ASparse;

    // vectors and arrays for particle generation
    std::vector<int> particlesInd;                               // vector holding total indices for each particle
    std::vector<int> particlesInterInd;                          // interfacial distance indices

    // initialize material properties for each node
    double density[SIZEA3];
    double heatCap[SIZEA3];
    double thermCond[SIZEA3];
    std::fill_n(density, SIZEA3, RHO0M);
    std::fill_n(heatCap, SIZEA3, KM);
    std::fill_n(thermCond, SIZEA3, CM);

    // initialize temperature
    std::vector< std::vector<double> > temperature;             // store temperature all time steps
    std::vector<double> theta(SIZEA3, THETA0);                  // store temperature for individual time steps

    // compute necessary matrices
    computeCoord(cubeCoord);                // compute the x-y-z coordinates
    computeBoundary(bNodes);                // find the boundary nodes
    computeAsparse(ASparse, bNodes);     // compute FDM mesh A matrix

    // compute laser profile
    double laserValues[SIZEA3];
    std::fill_n(laserValues, SIZEA3, 0.);
    std::cout << "test: " << cubeCoord[0][1] << std::endl;
    for (int i = 0; i < SIZEA3; i++){
        if (0.02 <= cubeCoord[i][1] <= 0.03 and 0.02 <= cubeCoord[i][2] <= 0.03){
            laserValues[i] = ABSORB * I0 * exp(-ABSORB * (LENGTH - cubeCoord[i][3]));
            std::cout << "laserValue: " << laserValues[i] << std::endl;
        }
    }

    // compute the random particles embedded in the material
    computeParticles(cubeCoord, particlesInd, particlesInterInd,
                     density, heatCap, thermCond);


    solutionScheme(ASparse, bNodes, temperature, theta,
                   density, heatCap, thermCond, laserValues);


//    write2file(cubeCoord, density); // write initial density result to output file


    return 0;
}

// Function definitions
void computeCoord(std::vector<std::vector< double>> &cubeCoord){
    /* structure cubeCoord:
     * cubeCoord = [node, x, y, z]
     */

    // required variables
    double coordDist = 0;               // distance incrementer
    std::vector<double> x;              // node numbering in x direction
    std::vector<double> y;              // node numbering in y direction
    std::vector<double> z;              // node numbering in z direction
    std::vector<double> xyz;            // node numbering in xyz direction

    std::cout << "--- Initializing computation --- " << std::endl;
    std::cout << "--- Constructing cubeCoord array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // populate vectors x-y-z
    for (int i = 0; i < NODE; i +=1 ){
        xyz.push_back(coordDist);
        coordDist += H;
    }

    // populate array coordCube with x-y-z directions
    int counter = 0;
    for (int i=0; i<NODE; i++){
        for (int j=0; j<NODE; j++){
            std::vector<double> temp;
            for (int k=0; k<NODE; k++){
                temp.push_back(counter);
                temp.push_back(xyz[k]);
                temp.push_back(xyz[j]);
                temp.push_back(xyz[i]);
//                std::cout << xyz[k] << " " << xyz[j]  << " " << xyz[i] << std::endl;
                counter++;
                cubeCoord.push_back(temp);
                temp.clear();
            }

        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "computational time: " << duration << "s\n" << std::endl;

}

void computeBoundary(std::vector< std::vector<int> >& bNodes){
    std::cout << "--- Initializing computation --- " << std::endl;
    std::cout << "--- Constructing bNodes array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<int> topBoundary;
    std::vector<int> topBoundary_h;
    std::vector<int> bottomBoundary;
    std::vector<int> bottomBoundary_h;
    std::vector<int> wallBoundary;

    int counter = 0;
    for (int i = 0; i < NODE; i++){
        for (int j = 0; j < NODE; j++){
            for (int k = 0; k < NODE; k++){
                if (i == 0){
                    // log the two boundary nodes
                    bottomBoundary.push_back(counter);              // boundary node
                    bottomBoundary_h.push_back(counter + SIZEA2);   // adjacent boundary node
                }
                if (i == NODE - 1){
                    topBoundary.push_back(counter);                 // boundary node
                    topBoundary_h.push_back(counter - SIZEA2);      // adjacent boundary node
                }
                if (j == 0 or j == NODE - 1 or k == 0 or k == NODE - 1){
                    wallBoundary.push_back(counter);
                }
                counter++;
            }
        }
    }
    bNodes.push_back(topBoundary);
    bNodes.push_back(topBoundary_h);
    bNodes.push_back(bottomBoundary);
    bNodes.push_back(bottomBoundary_h);
    bNodes.push_back(wallBoundary);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "compute time boundaries: " << duration << "s\n" << std::endl;
}

void computeAsparse(std::vector< std::vector<int> >& ASparse,
                    std::vector< std::vector<int> >& bNodes){
    /* Build a sparse A matrix
     *  A = [N^3, 4]
     *  [node, i, j, value]
     */
    std::cout << "--- Initializing computation --- " << std::endl;
    std::cout << "--- Constructing A sparse array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

//    int element = 0;
    for (int i = 0; i < SIZEA3; i++){
        if (SIZEA2 < i < (SIZEA3 - SIZEA2)){
            std::vector<int> tempVec;
            for (int j = 0; j < SIZEA3; j++){
                if ((j == i - SIZEA2) or (j == i - NODE) or (j == i - 1)){
//                    tempVec.push_back(element);
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(1);
                }else if (j == i){
//                    tempVec.push_back(element);
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(-6);
                }else if((j == i + 1) or (j == i + NODE) or (j == i + SIZEA2)){
//                    tempVec.push_back(element);
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(1);
                }else{
//                    element++;
                    continue;
                }
                ASparse.push_back(tempVec);
                tempVec.clear();
            }
        }
    }

    /*
     * enforce dirichlet boundary conditions along the walls
     * loop through boundary wall nodes bNodes[4][i] and set the value of
     * ASparse[node][3] equal to zero
     */
    for (int i = 0; i < bNodes[4].size(); i++){
        ASparse[bNodes[4][i]][3] = 0;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "compute time A matrix: " << duration << "s" << std::endl;
}

void computeParticles(std::vector< std::vector<double> >& cubeCoord,
                      std::vector<int> &particlesInd,
                      std::vector<int> &particlesInterInd,
                      double density[SIZEA3],
                      double heatCap[SIZEA3],
                      double thermCond[SIZEA3]){

    int nParticleNode = std::round(SIZEA3 * VP);     // total number of host nodes for particles
    double partDist, nodeParticle, randLoc, testVar;

    int counter1 = 0;
    while ((particlesInd.size() < nParticleNode) and (counter1 < 10000)){

        // choose a random node to generate particle
        // https://stackoverflow.com/questions/19665818/generate-random-numbers-using-c11-random-library
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::uniform_real_distribution<> dist{0, 1};
        randLoc = dist(gen);          // ensure diameter is positive

        // Generate random seed location for particle
        std::cout << "\n-- Generating particle " << counter1 + 1 << std::endl;
        nodeParticle = ceil(SIZEA3 * randLoc);

        // find nodes that are within the distance of the seed location
        for (int i = 0; i < SIZEA3; i++){

            // compute distance between generated particle and all other particles
            partDist = sqrt(pow(cubeCoord[i][1] - cubeCoord[nodeParticle][1], 2) +
                            pow(cubeCoord[i][2] - cubeCoord[nodeParticle][2], 2) +
                            pow(cubeCoord[i][3] - cubeCoord[nodeParticle][3], 2));

            // determine if distance is within the range of another particle
            if (partDist <= RPART){
                particlesInd.push_back(i);
            }else if (INTTHICK != 0 and partDist <= RPART + INTTHICK * H){
                particlesInterInd.push_back(i);

            }else{
                continue;
            }
        }

        uniqueVec(particlesInd);
        uniqueVec(particlesInterInd);

//        std::cout << "particleInd_.size(): " << particlesInd.size() << std::endl;
//        std::cout << "particleInterInd_.size(): " << particlesInterInd.size() << std::endl;

        // assign interfacial material properties
        for (int i = 0; i < particlesInterInd.size(); i++){
            density[particlesInterInd[i]] = (RHO0M + RHO0P) / 2.;
            heatCap[particlesInterInd[i]] = (CP + CM) / 2.;
            thermCond[particlesInterInd[i]] = (KP + KM) / 2.;
        }

        // assign particle material properties
        for (int i = 0; i < particlesInd.size(); i++){
            density[particlesInd[i]] = RHO0P;
            heatCap[particlesInd[i]] = CP;
            thermCond[particlesInd[i]] = KP;
        }
        counter1++;
    }
}


void solutionScheme(std::vector< std::vector<int> >& ASparse,
                    std::vector< std::vector<int> >& bNodes,
                    std::vector< std::vector<double> >& temperature,
                    std::vector<double>& theta,
                    double density[SIZEA3], double heatCap[SIZEA3],
                    double thermCond[SIZEA3], double laserValues[SIZEA3]){

    // begin time stepping
    double prefix1, prefix2;
    int row, col, val;
    int nnz = ASparse.size();
    for (int t = 0; t < SIZETIME; t++){
        std::cout << "time: " << t << std::endl;
        // coordinate-wise sparse-matrix-vector multiplication
        // http://www.mathcs.emory.edu/~cheung/Courses/561/Syllabus/3-C/sparse.html
        double AtimesTheta[SIZEA3];
        std::fill_n(AtimesTheta, SIZEA3, 0);
        for (int i = 0; i < nnz; i++){
            row = ASparse[i][0];
            col = ASparse[i][1];
            val = ASparse[i][2];
            AtimesTheta[row] = AtimesTheta[row] + val * theta[col];
        }

        // using solve for temperatures at the next time step
        for (int j = 0; j < SIZEA3; j++){
            prefix1 = DT / density[j] / heatCap[j];
            prefix2 = thermCond[j] / pow(H, 2.0);
            theta[j] = theta[j] + prefix1 * (prefix2 * AtimesTheta[j] + laserValues[j]);
        }

        // enforce Neumann boundary conditions on the top boundary
        for (int k = 0; k < bNodes[0].size(); k++){
            theta[bNodes[0][k]] = theta[bNodes[1][k]];
        }

        // enforce Neumann boundary conditions on the bottom boundary
        for (int l = 0; l < bNodes[3].size(); l++){
            theta[bNodes[3][l]] = theta[bNodes[4][l]];
        }

        // store temperature results
        if (t % 100){
            temperature.push_back(theta);
        }else{
            continue;
        }
    }

}

void write2file(std::vector< std::vector<double> > &cubeCoord,
                double density[SIZEA3]){

    // write to file
    printState.open("/Users/bhopro/Desktop/Berkeley/MSOL/Projects/Voxel/state.dat");
    printState << "X Y Z D" << std::endl;

    for (int i = 0; i < SIZEA3; i++){
        printState << cubeCoord[i][1] << " " << cubeCoord[i][2] << " " << cubeCoord[i][3]
                   << " " << density[i] << std::endl;
    }
    printState.close();
}

void uniqueVec(std::vector<int> &vec){
    // return only unique nodes in particleInd
    // https://www.geeksforgeeks.org/stdunique-in-cpp/

//    std::cout << "before: " << std::endl;
//    for (int i = 0; i < vec.size(); i++){
//        std::cout << vec[i] << " ";
//    }
//    std::cout << std::endl;

    // begin removing duplicate indices
    std::vector<int>::iterator ip;
    ip = std::unique(vec.begin(), vec.begin() + vec.size());
    vec.resize(std::distance(vec.begin(), ip));

//    std::cout << "after: " << std::endl;
//    for (ip = vec.begin(); ip != vec.end(); ++ip) {
//        std::cout << *ip << " ";
//    }
//    std::cout << std::endl;
}

void print2dVecDouble(std::vector< std::vector<double> > &vec){
    for (int i = 0; i < vec.size(); i++){
        for (int j = 0; j < vec[i].size(); j++){
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print2dVecInt(std::vector< std::vector<int> > &vec){
    for (int i=0; i<vec.size(); i++){
        for (int j=0; j<vec[i].size(); j++){
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}





