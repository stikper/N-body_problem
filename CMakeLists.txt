cmake_minimum_required(VERSION 3.26)
project(N-body_problem)

set(CMAKE_CXX_STANDARD 17)

add_executable(N-body_problem main.cpp
        Modules/MatVec.cpp
        Modules/MatVec.h
        Modules/DataOut.cpp
        Modules/DataOut.h
        Modules/Body.cpp
        Modules/Body.h
        Modules/Consts.h
        Modules/Compute.cpp
        Modules/Compute.h)
