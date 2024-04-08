cmake_minimum_required(VERSION 3.27)
project(Interaction_SAW_Models)

set(CMAKE_CXX_STANDARD 17)


add_executable(Interaction_SAW_Models
        include/Lattice.h
        Lattice.cpp
        main.cpp
        include/Model.h)


target_include_directories(Interaction_SAW_Models PUBLIC include)
