cmake_minimum_required(VERSION 2.8.9)

# Use Eigen3.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)

# Base compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O3 -march=native -Wall") # Optimize for this machine, and turn on compiler warnings.
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")             # Optimize for this architecture.
endif()

# Enable PIC for shared library generation.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Include directories.
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(.)

add_executable(rvat rvat.cpp)
target_link_libraries(rvat vortexje)
