/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel3.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOMPILEVEL3_H_
#define IOMPILEVEL3_H_

#include "HeatTransfer.h"
#include "Settings.h"
#include "FileView.h"

#include <vector>
#include <mpi.h>

class IOmpiLevel3
{
 public:
  IOmpiLevel3() = default;
  
  IOmpiLevel3( const Settings& s, MPI_Comm communicator );
  
  ~IOmpiLevel3() = default;
  
  IOmpiLevel3( IOmpiLevel3 const& other ) = delete;
  
  IOmpiLevel3( IOmpiLevel3&& other ) = default;
  
  IOmpiLevel3& operator=( IOmpiLevel3 const& other ) = delete;
  
  IOmpiLevel3& operator=( IOmpiLevel3&& other ) noexcept {
    IOmpiLevel3 tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  }
  
  void swap( IOmpiLevel3& other ) noexcept {
    using std::swap;
    swap( _filehandle, other._filehandle );
    swap( _fileview, other._fileview );
    swap( _communicator, other._communicator );
    swap( _outputfilename, other._outputfilename );
    swap( _buffercount, other._buffercount );
    swap( _rank, other._rank );
    swap( _nprocs, other._nprocs );
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
  MPI_File _filehandle;
  FileView _fileview;
  MPI_Comm _communicator;
  std::string _outputfilename;
  int _buffercount;
  int _rank;
  int _nprocs;
};

void swap( IOmpiLevel3& a, IOmpiLevel3& b ) noexcept;

#endif /* IOMPILEVEL3_H_ */
