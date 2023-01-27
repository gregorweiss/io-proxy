/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel0.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 *
 *  This refers to the example of Figure 7.2 in 'Using Advance MPI' of Gropp et al.
 */

#include "IOmpiLevel0.h"
#include "helper.h"

#include <cstdio>

#include <mpi.h>

IOmpiLevel0::IOmpiLevel0( const Settings& s, MPI_Comm comm )
  : _communicator{ comm }
  , _outputfilename{ MakeFilename( s.outputfile, ".mpiio_write" ) }
  , _buffercount{ static_cast<int>( s.ndx * s.ndy ) }
  , _rank{ getRank( comm ) }
  , _nprocs{ getNProcs( comm ) } {}

void IOmpiLevel0::write( int step,
                         const HeatTransfer& ht,
                         const Settings& s,
                         MPI_Comm comm ) {
  std::vector<double> v = ht.data_noghost();

  _outputfilename = MakeFilename( s.outputfile, ".mpiio_write", -1, step );
  
  // Open file and set initial rank-related offset
  MPI_File_open( comm,
                 _outputfilename.c_str(),
                 MPI_MODE_CREATE | MPI_MODE_WRONLY,
                 MPI_INFO_NULL,
                 &_filehandle );
  MPI_Offset offset = _rank * _buffercount * sizeof( double );
  MPI_File_seek( _filehandle, offset, MPI_SEEK_SET );
  
  MPI_File_write( _filehandle,
                  v.data(),
                  _buffercount,
                  MPI_DOUBLE,
                  MPI_STATUS_IGNORE);
  
  MPI_File_close( &_filehandle );
}

void IOmpiLevel0::read( const int step,
                        std::vector<double>& buffer,
                        const Settings& s,
                        MPI_Comm comm ) {
  _outputfilename = MakeFilename( s.outputfile, ".mpiio_write", -1, step );

  // Open file and set initial rank-related offset
  MPI_File_open( comm,
                 _outputfilename.c_str(),
                 MPI_MODE_RDONLY,
                 MPI_INFO_NULL,
                 &_filehandle );
  MPI_Offset offset = _rank * _buffercount * sizeof( double );
  MPI_File_seek( _filehandle, offset, MPI_SEEK_SET );
  
  MPI_File_read( _filehandle, buffer.data(), _buffercount, MPI_DOUBLE, MPI_STATUS_IGNORE);
  
  MPI_File_close( &_filehandle );
}

void IOmpiLevel0::remove( const int step ) {
  if ( _rank == 0 )
    std::remove( _outputfilename.c_str());
}

void swap( IOmpiLevel0& a, IOmpiLevel0& b ) noexcept {
  a.swap( b );
}
