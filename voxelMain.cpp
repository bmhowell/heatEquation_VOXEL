//
// Voxel code CPP
// Created by Brian Howell on 8/17/21.
//


#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <chrono>
#include <random>
#include <vector>

// particle parameters
const double THETAW = 300.;                             // |    K    |  temperature at the wall
const double THETA0 = 300.;                             // |    K    |  initial temperature
const double I0 = 2 * pow(10, 8);         // |  W/m^2  |  initial laser intensity
const double ABSORB = 25.0;                             // |   1/m   |  absorption constant
const double LENGTH = 0.05;                             // |    m    |  sample length
const int NODE = 21;                                    // |   ___   |  number of nodes

// material parameters
const int KP = 10;                                      // |  W/m-K  |  thermal conductivity
const int C = 450;                                      // | J/kg-K  |  heat capacity
const int RHO0 = 3000;                                  // | kg/m^3  |  initial material density
const double VP = 0.3;                                  // |   ---   |  volume fraction of particles
const double RPART = LENGTH / 10.;                      // |    m    |  radius of the particles

// simulation parameters
const double TF = 2.;                                   // |    s    |  final simulation time
const double DT = pow(10, -4);            // |    s    |  initial time step

const int SIZEA3 = (int) pow(NODE, 3);
const int SIZEA2 = (int) pow(NODE, 2);

// data outputs
std::ofstream printParticles;
std::ofstream printDiameters;

// function declarations
void computeCoord(std::vector< std::vector< double > >&);
// computeCoord - compute the coords in the x-y-z directions
// @param vector< vector<double> > - the 2D vector to be filled
// @param double N -  Number of nodes
// @param double L -  Sample length

void computeBoundary(std::vector< std::vector<double> >&);
// boundaryNodes - computes which nodes are on the sides, top and bottom of the cube
// @params vector< vector<double> > - input 2D vector to be filled
// @params double - number of nodes
// @return vector< vector<double> > - vector of vectors containing each side

void computeAsparse(std::vector< std::vector<int> >&);

void print2dVecInt(std::vector< std::vector<int> >&);
// print2dVec - prints 2D vector

int main(){
    std::cout << SIZEA3 << std::endl;
    // initialize cubeCoord, bNode, and A matrix vectors
    std::vector< std::vector<double> > cubeCoord;
    std::vector< std::vector<double> > bNodes;
    std::vector< std::vector<int> > ASparse;

    computeCoord(cubeCoord);
    computeBoundary(bNodes);
    computeAsparse(ASparse);


    return 0;
}

// Function definitions
void computeCoord(std::vector<std::vector< double>> &cubeCoord){

    // required variables
    double h = LENGTH / (NODE - 1);     // physical discretization length
    double coordDist = 0;               // distance incrementer
    std::vector<double> x;              // node numbering in x direction
    std::vector<double> y;              // node numbering in y direction
    std::vector<double> z;              // node numbering in z direction
    std::vector<double> xyz;            // node numbering in xyz direction

    std::cout << "--- Initializing computation --- " << std::endl;
    std::cout << "--- Constructing cubeCoord array ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // populate vectors x-y-z
    for (int i=0; i<NODE; i+=1){
        xyz.push_back(coordDist);
        coordDist += h;
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
    std::cout << "computational time: " << duration << "s" << std::endl;

}

void computeBoundary(std::vector< std::vector<double> >& bNodes){

    std::vector<double> topBoundary;
    std::vector<double> topBoundary_h;
    std::vector<double> bottomBoundary;
    std::vector<double> bottomBoundary_h;
    std::vector<double> wallBoundary;

    int counter = 0;
    for (int i=0; i<NODE; i++){
        for (int j=0; j<NODE; j++){
            for (int k=0; k<NODE; k++){
                if (i==0){
                    bottomBoundary.push_back(counter);
                    bottomBoundary_h.push_back(counter + pow(NODE, 2.0));
                }
                if (i==NODE-1){
                    topBoundary.push_back(counter);
                    topBoundary_h.push_back(counter);
                }
                if (j==0 or j==NODE-1 or k==0 or k==NODE-1){
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
}

void computeAsparse(std::vector< std::vector<int> >& ASparse){
    /* Build a sparse A matrix
     *  A = [N^3, 4]
     *  [node, i, j, value]
     */
    auto start = std::chrono::high_resolution_clock::now();

    int node = 0;
    for (int i=0; i<SIZEA3; i++){
        if (SIZEA2 < i < (SIZEA3 - SIZEA2)){
            std::vector<int> tempVec;
            for (int j=0; j<SIZEA3; j++){
                if ((j == i - SIZEA2) or (j == i - NODE) or (j == i - 1)){
                    tempVec.push_back(node);
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(1);
                }else if (j==i){
                    tempVec.push_back(node);
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(-6);
                }else if((j == i + 1) or (j == i + NODE) or (j == i + SIZEA2)){
                    tempVec.push_back(node);
                    tempVec.push_back(i);
                    tempVec.push_back(j);
                    tempVec.push_back(1);
                }else{
                    node++;
                    continue;
                }
                ASparse.push_back(tempVec);
                tempVec.clear();
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(stop - start)).count() / 1e6;
    std::cout << "computational time: " << duration << "s" << std::endl;
}

void print2dVecInt(std::vector< std::vector<int> > &vec){
    for (int i=0; i<vec.size(); i++){
        for (int j=0; j<vec[i].size(); j++){
            std::cout << vec[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print2dVecInt(std::vector< std::vector<int> > &vec);


