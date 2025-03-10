/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel0.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOMPILEVEL0_H_
#define IOMPILEVEL0_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <vector>
#include <mpi.h>

class IOmpiLevel0
{
 public:
  IOmpiLevel0() = default;
  
  IOmpiLevel0( const Settings& s, MPI_Comm comm );
  
  ~IOmpiLevel0() = default;
  
  IOmpiLevel0( IOmpiLevel0 const& other ) = delete;
  
  IOmpiLevel0& operator=( IOmpiLevel0 const& other ) = delete;
  
  IOmpiLevel0( IOmpiLevel0&& other ) = default;
  
  IOmpiLevel0& operator=( IOmpiLevel0&& other ) noexcept {
    IOmpiLevel0 tmp{ std::move( other ) };
    swap( other );
    return *this;
  }
  
  void swap( IOmpiLevel0& other ) noexcept {
    using std::swap;
    swap( _filehandle, other._filehandle );
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
  MPI_File _filehandle{};
  MPI_Comm _communicator{};
  std::string _outputfilename{};
  MPI_Offset _disp{};
  int _buffercount{};
  int _rank{};
  int _nprocs{};
};

void swap( IOmpiLevel0& a, IOmpiLevel0& b ) noexcept;

#endif /* IOMPILEVEL0_H_ */
