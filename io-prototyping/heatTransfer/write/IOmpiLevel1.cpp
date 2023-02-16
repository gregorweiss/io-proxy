/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel1.cpp
 *
 *  Created on: Feb 2023
 *      Author: Gregor Weiss
 *
 *  This refers to the example of Figure 7.2 in 'Using Advance MPI' of Gropp et al.
 */

#include "IOmpiLevel1.h"
#include "helper.h"

#include <cstdio>

#include <mpi.h>

IOmpiLevel1::IOmpiLevel1( const Settings& s, MPI_Comm comm )
  : _communicator{ comm }
  , _outputfilename{ MakeFilename( s.outputfile, ".mpiio_write" ) }
  , _disp{ static_cast<MPI_Offset>( sizeof(double) * s.gndx * s.gndy ) }
  , _buffercount{ static_cast<int>( s.ndx * s.ndy ) }
  , _rank{ getRank( comm ) }
  , _nprocs{ getNProcs( comm ) } {}

void IOmpiLevel1::write( int step,
                         const HeatTransfer& ht,
                         const Settings& s,
                         MPI_Comm comm ) {
  _outputfilename = MakeFilename( s.outputfile, ".mpiio_write", -1, step );
  
  // Open file and set initial rank-related offset
  MPI_File_open( comm,
                 _outputfilename.c_str(),
                 MPI_MODE_CREATE | MPI_MODE_WRONLY,
                 MPI_INFO_NULL,
                 &_filehandle );

  MPI_Offset offset = _rank * _buffercount * sizeof( double );
  for ( const auto& iteration : ht.m_TIterations ) {
    MPI_File_seek( _filehandle, offset, MPI_SEEK_SET );
    MPI_File_write_all( _filehandle,
                        iteration.data(),
                        _buffercount,
                        MPI_DOUBLE,
                        MPI_STATUS_IGNORE);
    offset += _disp;
  }
  
  MPI_File_close( &_filehandle );
}

void IOmpiLevel1::read( const int step,
                        std::vector<std::vector<double> >& buffer,
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
  for ( auto& iteration : buffer ) {
    MPI_File_seek( _filehandle, offset, MPI_SEEK_SET );
    MPI_File_read_all( _filehandle,
                       iteration.data(),
                       _buffercount,
                       MPI_DOUBLE,
                       MPI_STATUS_IGNORE);
    offset += _disp;
  }
  
  MPI_File_close( &_filehandle );
}

void IOmpiLevel1::remove( const int step ) {
  if ( _rank == 0 )
    std::remove( _outputfilename.c_str());
}

void swap( IOmpiLevel1& a, IOmpiLevel1& b ) noexcept {
  a.swap( b );
}
