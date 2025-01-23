cmake_minimum_required(VERSION 3.10)
project(finite_element LANGUAGES CXX CUDA)

# -------------------------------------------------------------------
# 1. Set C++ and CUDA Standards
# -------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# -------------------------------------------------------------------
# 2. Define Build Options
# -------------------------------------------------------------------
option(USE_CUDA "Enable CUDA support" ON)

# -------------------------------------------------------------------
# 3. Locate Kokkos
# -------------------------------------------------------------------
# Set Kokkos_DIR to the directory containing KokkosConfig.cmake
# Adjust the path if Kokkos is installed elsewhere
set(Kokkos_DIR "${HOME}/kokkos/install/lib/cmake/Kokkos")

find_package(Kokkos REQUIRED PATHS "${Kokkos_DIR}" NO_DEFAULT_PATH)

# -------------------------------------------------------------------
# 4. Set Compilation Flags
# -------------------------------------------------------------------
# Equivalent to CXXFLAGS in Makefile
set(CXX_FLAGS "-O3 -DREGIME_3D -DSTARTHALF -use_fast_math -I./include/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS}")

# If using CUDA, set additional CUDA flags
if(USE_CUDA)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -DREGIME_3D -DSTARTHALF -use_fast_math")
endif()

# -------------------------------------------------------------------
# 5. Gather Source Files
# -------------------------------------------------------------------
file(GLOB SRC_FILES "*.cpp")

# -------------------------------------------------------------------
# 6. Add Executable
# -------------------------------------------------------------------
add_executable(${PROJECT_NAME} ${SRC_FILES})

# -------------------------------------------------------------------
# 7. Link Kokkos
# -------------------------------------------------------------------
target_link_libraries(${PROJECT_NAME} Kokkos::kokkos)

# -------------------------------------------------------------------
# 8. Optional: Define Macros or Additional Includes
# -------------------------------------------------------------------
# If you have additional include directories or libraries, add them here
# target_include_directories(${PROJECT_NAME} PRIVATE /path/to/includes)
# target_link_libraries(${PROJECT_NAME} PRIVATE /path/to/libs)

# -------------------------------------------------------------------
# 9. Set CUDA Properties (Optional)
# -------------------------------------------------------------------
if(USE_CUDA)
    set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # Add any CUDA-specific compile options or properties here
endif()
