#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(IO-Prototyping C CXX)

# add cmake modules: for all `include(...)` first look here
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED)

# Set the ADIOS2 library name
set(ADIOS2_LIB adios2)
if (ADIOS2_VERSION_MAJOR GREATER_EQUAL 2)
  if (ADIOS2_VERSION_MINOR GREATER_EQUAL 6)
    set(ADIOS2_LIB adios2::cxx11_mpi)
  endif ()
endif ()

# These includes publish function names.
include(ProcessOptions)

# Link with SIONlib
process_with_sionlib()
foreach (def ${SIONLIB_DEFINES})
    set(ALL_CFLAGS "${ALL_CFLAGS} ${def}")
    set(ALL_CXXFLAGS "${ALL_CXXFLAGS} ${def}")
endforeach ()

set(ALL_LIBS "${SIONLIB_LIBRARIES}")
set(ALL_INCLUDES "${ALL_INCLUDES} ${SIONLIB_INCLUDE}")

add_subdirectory(write)

