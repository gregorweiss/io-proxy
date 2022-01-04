#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(IO-Prototyping C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED)

# Set the ADIOS2 library name
set(ADIOS2_LIB adios2)
if(ADIOS2_VERSION GREATER_EQUAL 2.6)
    set(ADIOS2_LIB adios2::cxx11_mpi)
endif()

add_subdirectory(write)
add_subdirectory(read)
#add_subdirectory(read_fileonly)
#add_subdirectory(inline)
