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

#include <mpi.h>

FileView::FileView( const Settings& settings, MPI_Comm communicator )
  : _communicator{ communicator }
  , _psizes{ static_cast<int>(settings.npx), static_cast<int>(settings.npy) }
  , _gsizes{ static_cast<const int>(settings.gndx), static_cast<const int>(settings.gndy) }
  , _globsizes{ static_cast<int>(settings.gndx * settings.gndy) }
  , _rank{ getRank( _communicator ) }
  , _nprocs{ getNProcs( _communicator ) } {
  initialize();
}

FileView::FileView( FileView const& other )
  : _communicator{ other._communicator }
  , _psizes{ other._psizes[0], other._psizes[1] }
  , _gsizes{ other._gsizes[0], other._gsizes[1] }
  , _globsizes{ other._globsizes }
  , _rank{ other._rank }
  , _nprocs{ other._nprocs } {
  initialize();
}

FileView::~FileView() noexcept {
  if ( _filetype )
    MPI_Type_free( &_filetype );
  if ( _darraytype )
    MPI_Type_free( &_darraytype );
}

void FileView::initialize() {
  MPI_Type_create_darray( _nprocs,
                          _rank,
                          2,
                          _gsizes,
                          _distribs,
                          _dargs,
                          _psizes,
                          MPI_ORDER_C,
                          MPI_DOUBLE,
                          &_darraytype );
  MPI_Aint lowerbound{ 0 };
  MPI_Aint extend = static_cast<MPI_Aint>( _globsizes * sizeof( double ));
  MPI_Type_create_resized( _darraytype, lowerbound, extend, &_filetype );
  MPI_Type_commit( &_filetype );
}

void swap( FileView& a, FileView& b ) noexcept {
  a.swap( b );
}

IOmpiLevel3::IOmpiLevel3( const Settings& s, MPI_Comm communicator )
  : _fileview{ s, communicator }
  , _communicator{ communicator }
  , _outputfilename{ MakeFilename( s.outputfile, "mpi_write_all" ) }
  , _buffercount{ static_cast<int>( s.ndx * s.ndy ) }
  , _rank{ getRank( _communicator ) }
  , _nprocs{ getNProcs( _communicator ) } {}

void IOmpiLevel3::write( int step,
                         const HeatTransfer& ht,
                         const Settings& s,
                         MPI_Comm comm ) {
  std::vector<double> v = ht.data_noghost();

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
                     MPI_INFO_NULL);
  MPI_File_write_all( filehandle_onestep,
                      v.data(),
                      _buffercount,
                      MPI_DOUBLE,
                      MPI_STATUS_IGNORE);
  
  MPI_File_close( &filehandle_onestep );
}


void IOmpiLevel3::read( const int step,
                        std::vector<double>& buffer,
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
                     MPI_INFO_NULL);
  MPI_File_read_all( filehandle_onestep,
                     buffer.data(),
                     _buffercount,
                     MPI_DOUBLE,
                     MPI_STATUS_IGNORE);

  MPI_File_close( &filehandle_onestep );
}

void IOmpiLevel3::remove( const int step ) {
  if ( _rank == 0 )
    std::remove( _outputfilename.c_str());
}

void swap( IOmpiLevel3& a, IOmpiLevel3& b ) noexcept {
  a.swap( b );
}
