#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_executable(heatTransfer_write_adios2
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_adios2.cpp
)
target_link_libraries(heatTransfer_write_adios2 
  ${ADIOS2_LIB} MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})

add_executable(heatTransfer_write_ascii
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_ascii.cpp
)
target_link_libraries(heatTransfer_write_ascii
  MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
)

add_executable(heatTransfer_write_binary
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_binary.cpp
)
target_link_libraries(heatTransfer_write_binary
  MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
)

add_executable(heatTransfer_write_mpiio_write
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_mpiio_write.cpp
)
target_link_libraries(heatTransfer_write_mpiio_write
  MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
)

add_executable(heatTransfer_write_mpiio_write_at
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_mpiio_write_at.cpp
)
target_link_libraries(heatTransfer_write_mpiio_write_at
  MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
)

add_executable(heatTransfer_write_mpiio_write_all
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_mpiio_write_all.cpp
)
target_link_libraries(heatTransfer_write_mpiio_write_all
  MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
)

add_executable(heatTransfer_write_mpiio_write_all_indexed
  main.cpp
  HeatTransfer.cpp
  Settings.cpp
  IO_mpiio_write_all_indexed.cpp
)
target_link_libraries(heatTransfer_write_mpiio_write_all_indexed
  MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
)

if(ADIOS2_HAVE_HDF5)
  add_executable(heatTransfer_write_hdf5
    main.cpp
    HeatTransfer.cpp
    Settings.cpp
    IO_hdf5_a.cpp
  )
  target_include_directories(heatTransfer_write_hdf5
    PRIVATE ${HDF5_C_INCLUDE_DIRS}
  )
  target_link_libraries(heatTransfer_write_hdf5
    ${HDF5_C_LIBRARIES} MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
  )

  if(HDF5_IS_PARALLEL)
    add_executable(heatTransfer_write_ph5
      main.cpp
      HeatTransfer.cpp
      Settings.cpp
      IO_ph5.cpp
    )
    target_include_directories(heatTransfer_write_ph5
      PRIVATE ${HDF5_C_INCLUDE_DIRS}
    )
    target_link_libraries(heatTransfer_write_ph5
      ${HDF5_C_LIBRARIES} MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
    )
  endif()

  if(NOT (HDF5_VERSION VERSION_LESS 1.11))
    add_executable(heatTransfer_write_h5mixer
      main.cpp
      HeatTransfer.cpp
      Settings.cpp
      IO_h5mixer.cpp
    ) 
    target_link_libraries(heatTransfer_write_h5mixer
      ${ADIOS2_LIB} MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT}
    )
  endif()
endif()
