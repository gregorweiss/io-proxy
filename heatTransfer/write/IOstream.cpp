/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOstream.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#include "IOstream.h"
#include "helper.h"

#include <iostream> //std::cout

#include <stdio.h>
#include <unistd.h>

#include <cassert>

IOstream::IOstream( const Settings& s, MPI_Comm comm )
  : _filename{ MakeFilename( s.outputfile, ".dat", s.rank ) } {}

IOstream::~IOstream() {}

void IOstream::write( int step,
                      const HeatTransfer& ht,
                      const Settings& s,
                      MPI_Comm comm ) {
  _filename = MakeFilename( s.outputfile, ".dat", s.rank, step );
  _filestream = fopen( _filename.c_str(), "w" );
  
  auto write_size = s.ndx * s.ndy;
  for ( const auto& iteration : ht.m_TIterations ) {
    fwrite( reinterpret_cast<const char*>(iteration.data()),
            sizeof( double ),
            write_size,
            _filestream );
  }
  
  fclose( _filestream );
  fsync( fileno( _filestream ));
}

void IOstream::read( const int step,
                     std::vector<std::vector<double> >& buffer,
                     const Settings& s,
                     MPI_Comm comm ) {
  _filename = MakeFilename( s.outputfile, ".dat", s.rank, step );
  _filestream = fopen( _filename.c_str(), "r" );

  auto read_size = s.ndx * s.ndy;
  for ( auto& iteration : buffer ) {
    size_t count = fread( reinterpret_cast<char*>(iteration.data()),
                          sizeof( double ),
                          read_size,
                          _filestream );
    assert(count==read_size);
  }

  fclose( _filestream );
}

void IOstream::remove( const int step ) {
  std::remove( _filename.c_str());
}
