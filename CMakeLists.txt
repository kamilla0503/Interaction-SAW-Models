cmake_minimum_required(VERSION 3.27)
project(Interaction_SAW_Models)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-DREGIME_2D -DCHEKMAP2D)
add_definitions(-DSTARTDEFAULT)
add_definitions(-DSEED)
#add_definitions(-DXYSHORT)
add_executable(Interaction_SAW_Models
        #include/Lattice.h
        Model.cpp
        Lattice.cpp
        Monte_Carlo.cpp
        main.cpp
        include/ExactEnumeration.h
        ExactEnumeration.cpp
        #include/Model.h
        #include/MonteCarlo.h
)


target_include_directories(Interaction_SAW_Models PUBLIC include)
