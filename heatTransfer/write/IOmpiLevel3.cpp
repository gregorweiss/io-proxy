/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel3.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 *
 *  This refers to the example of Figure 7.2 in 'Using Advance MPI' of Gropp et al.
 */

#include "IOmpiLevel3.h"
#include "helper.h"

#include <cstdio>
#include <iostream> // cout

#include <mpi.h>

auto choose_layout( std::string_view format ) {
  std::function<void(const Settings&, MPI_Datatype&)> gen_filetype{ nullptr };
  if ( format.find("1Dsubarray") != std::string::npos ) {
    gen_filetype = subarray1D;
  } else if ( format.find("2Dsubarray") != std::string::npos ) {
    if ( format.find("contiguous") != std::string::npos ) {
      gen_filetype = subarray2D_contiguous;
    } else {
      gen_filetype = subarray2D;
    }
  } else if ( format.find("1Ddarray") != std::string::npos ) {
    gen_filetype = darray1D;
  } else if ( format.find("2Ddarray") != std::string::npos ) {
    if ( format.find("contiguous") != std::string::npos ) {
      gen_filetype = darray2D_contiguous;
    } else {
      gen_filetype = darray2D;
    }
  } else {
    throw std::invalid_argument("Choose a file view form: 1Dsubarray, 1Ddarray, 2Ddarray.");
  }
  return gen_filetype;
}

IOmpiLevel3::IOmpiLevel3( const Settings& s, MPI_Comm communicator )
  : _fileview{ s, communicator, choose_layout( s.format ) }
  , _communicator{ communicator }
  , _outputfilename{ MakeFilename( s.outputfile, "mpi_write_all" ) }
  , _buffercount{ static_cast<int>( s.ndx * s.ndy * s.ndz ) }
  , _rank{ getRank( _communicator ) }
  , _nprocs{ getNProcs( _communicator ) } {}

void IOmpiLevel3::write( int step,
                         const HeatTransfer& ht,
                         const Settings& s,
                         MPI_Comm comm ) {
  _outputfilename = MakeFilename( s.outputfile, ".mpiio_write_all", -1, step );
  
  // Open file and set file view
  MPI_File filehandle_onestep;
  MPI_File_open( comm,
                 _outputfilename.c_str(),
                 MPI_MODE_CREATE | MPI_MODE_WRONLY,
                 MPI_INFO_NULL,
                 &filehandle_onestep );

  MPI_File_set_view( filehandle_onestep,
                     0,
                     MPI_DOUBLE,
                     _fileview._filetype,
                     "native",
                     MPI_INFO_NULL );
  for ( const auto& iteration : ht.m_TIterations ) {
    MPI_File_write_all( filehandle_onestep,
                        iteration.data(),
                        _buffercount,
                        MPI_DOUBLE,
                        MPI_STATUS_IGNORE );
  }
  
  MPI_File_close( &filehandle_onestep );
}

void IOmpiLevel3::read( const int step,
                        std::vector<std::vector<double> >& buffer,
                        const Settings& s,
                        MPI_Comm comm ) {
  _outputfilename = MakeFilename( s.outputfile, ".mpiio_write_all", -1, step );

  // Open file and set file view
  MPI_File filehandle_onestep;
  MPI_File_open( comm,
                 _outputfilename.c_str(),
                 MPI_MODE_RDONLY,
                 MPI_INFO_NULL,
                 &filehandle_onestep );

  MPI_File_set_view( filehandle_onestep,
                     0,
                     MPI_DOUBLE,
                     _fileview._filetype,
                     "native",
                     MPI_INFO_NULL );
  for ( auto& iteration : buffer ) {
    MPI_File_read_all( filehandle_onestep,
                       iteration.data(),
                       _buffercount,
                       MPI_DOUBLE,
                       MPI_STATUS_IGNORE );
  }

  MPI_File_close( &filehandle_onestep );
}

void IOmpiLevel3::remove( const int step ) {
  if ( _rank == 0 )
    std::remove( _outputfilename.c_str());
}

void swap( IOmpiLevel3& a, IOmpiLevel3& b ) noexcept {
  a.swap( b );
}
