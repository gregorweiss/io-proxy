/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOsion.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifdef HAVE_SIONLIB

#include "IOsion.h"
#include "helper.h"

#include <stdio.h>

IOsion::IOsion( const Settings& s, MPI_Comm comm )
  : _fileName{ MakeFilename( s.outputfile, ".dat", s.rank ) }
  , _communicator{ comm }
  , _chunkSize{ 10 * 1024 * 1024 }
  , _fsBlockSize{ -1 }
  , _numFiles{ 1 }
  , _rank{ getRank( _communicator ) }
  , _sionFileId{ 0 }
  , _filePtr{ nullptr }
  , _newFileName{ nullptr } {}

void IOsion::write( int step,
                    const HeatTransfer& ht,
                    const Settings& s,
                    MPI_Comm comm ) {
  std::vector<double> v = ht.data_noghost();

  _fileName = MakeFilename( s.outputfile, ".sion", -1, step );
  _sionFileId = sion_paropen_mpi( _fileName.c_str(), "bw", &_numFiles, _communicator, &_communicator,
                                   &_chunkSize, &_fsBlockSize, &_rank, &_filePtr, &_newFileName );
  
  sion_fwrite( v.data(),
               sizeof( double ),
               s.ndx * s.ndy,
               _sionFileId );
  
  sion_parclose_mpi( _sionFileId );
}

void IOsion::read( const int step,
                   std::vector<double>& ht,
                   const Settings& s,
                   MPI_Comm comm ) {
    std::cout << "IOsion::read not implemented for sion format." << std::endl;
}

void IOsion::remove( const int step ) {
  std::remove( _fileName.c_str());
}

#endif
