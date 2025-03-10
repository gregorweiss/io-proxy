#
#  Distributed under the OSI-approved Apache License, Version 2.0.  See
#  accompanying file Copyright.txt for details.
#
#  thirdparty.cmake
#
#  Created on: 2025
#      Author: Gregor Weiss
#
#

set(ADIOS2_CHECKOUT_VERSION
    "2.10.2"
    CACHE STRING "Use specific version of ADIOS2")
mark_as_advanced(ADIOS2_CHECKOUT_VERSION)

find_package(ADIOS2 ${ADIOS2_CHECKOUT_VERSION} QUIET)

include(cmake/CPM.cmake)

if(NOT ${ADIOS2_FOUND})
  set(ADIOS2_OPTIONS
      "BUILD_TYPE Release"
      "ADIOS2_USE_Fortran OFF"
      "ADIOS2_USE_Python OFF"
      "ADIOS2_BUILD_EXAMPLES OFF"
      "BUILD_TESTING OFF"
      "ADIOS2_USE_Profiling OFF")

  cpmaddpackage(
    NAME
    adios2
    GITHUB_REPOSITORY
    ornladios/ADIOS2
    VERSION
    2.10.2
    OPTIONS
    ${ADIOS2_OPTIONS}
    ${ADIOS2_CUDA_OPTIONS}
    SYSTEM)
else()
  message(STATUS "Found ADIOS2 ${ADIOS2_CHECKOUT_VERSION}")
endif()
