#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_executable(heatTransfer_write
        main.cpp
        HeatTransfer.cpp
        Settings.cpp
        helper.cpp
        IOascii.cpp
        IOadios2.cpp
        IObinary.cpp
        IOmpiLevel0.cpp
        IOmpiLevel1.cpp
        IOmpiLevel3.cpp
        IOsion.cpp
        IOstream.cpp
        )
target_link_libraries(heatTransfer_write
        ${ADIOS2_LIB}
        MPI::MPI_C
        ${CMAKE_THREAD_LIBS_INIT}
        ${SIONLIB_LIBRARIES})
