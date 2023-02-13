/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOsion.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOSION_H_
#define IOSION_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <fstream>
#include <vector>
#include <mpi.h>

#include "sion.h"

class IOsion
{
 public:
  IOsion() = default;
  
  IOsion( const Settings& s, MPI_Comm comm );
  
  ~IOsion() = default;
  
  IOsion( IOsion const& other ) = delete;
  
  IOsion& operator=( IOsion const& other ) = delete;
  
  IOsion( IOsion&& other ) = default;
  
  IOsion& operator=( IOsion&& other ) noexcept {
    IOsion tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  }
  
  void swap( IOsion& other ) noexcept {
    using std::swap;
    swap( _fileName, other._fileName );
    swap( _communicator, other._communicator );
    swap( _chunkSize, other._chunkSize );
    swap( _fsBlockSize, other._fsBlockSize );
    swap( _rank, other._rank );
    swap( _sionFileId, other._sionFileId );
    swap( _filePtr, other._filePtr );
    swap( _newFileName, other._newFileName );
  }
  
  void write( int step,
              const HeatTransfer& ht,
              const Settings& s,
              MPI_Comm comm );
  
  void read( const int step,
             std::vector<std::vector<double> >& buffer,
             const Settings& s,
             MPI_Comm comm );

  void remove( const int step );
 
 private:
  std::string _fileName{};
  MPI_Comm _communicator;
  sion_int64 _chunkSize{ 10 * 1024 * 1024 };
  sion_int32 _fsBlockSize{ -1 };
  int _numFiles{ 1 };
  int _rank{};
  int _sionFileId{};
  FILE* _filePtr{ nullptr };
  char* _newFileName{ nullptr };
};

#endif /* IOSION_H_ */
