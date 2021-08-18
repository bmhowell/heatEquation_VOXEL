cmake_minimum_required(VERSION 3.19)
project(Voxel)

set(CMAKE_CXX_STANDARD 14)

add_executable(Voxel
        voxelMain.cpp helpers.h helpers.cpp)
