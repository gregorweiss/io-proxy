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

IOstream::IOstream( const Settings& s, MPI_Comm comm )
  : _filename{ MakeFilename( s.outputfile, ".dat", s.rank ) } {}

IOstream::~IOstream() {}

void IOstream::write( int step,
                      const HeatTransfer& ht,
                      const Settings& s,
                      MPI_Comm comm ) {
  std::vector<double> v = ht.data_noghost();

  _filename = MakeFilename( s.outputfile, ".dat", s.rank, step );
  _filestream = fopen( _filename.c_str(), "w" );
  
  fwrite( reinterpret_cast<const char*>(v.data()),
          sizeof( double ),
          s.ndx * s.ndy,
          _filestream );
  
  fclose( _filestream );
  fsync( fileno( _filestream ));
}

void IOstream::read( const int step,
                   std::vector<double>& ht,
                   const Settings& s,
                   MPI_Comm comm ) {
    std::cout << "IOstream::read not implemented for stream format." << std::endl;
}

void IOstream::remove( const int step ) {
  std::remove( _filename.c_str());
}
