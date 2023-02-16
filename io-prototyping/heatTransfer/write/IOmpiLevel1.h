/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel1.h
 *
 *  Created on: Feb 2023
 *      Author: Gregor Weiss
 */

#ifndef IOMPILEVEL1_H_
#define IOMPILEVEL1_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <vector>
#include <mpi.h>

class IOmpiLevel1
{
 public:
  IOmpiLevel1() = default;
  
  IOmpiLevel1( const Settings& s, MPI_Comm comm );
  
  ~IOmpiLevel1() = default;
  
  IOmpiLevel1( IOmpiLevel1 const& other ) = delete;
  
  IOmpiLevel1& operator=( IOmpiLevel1 const& other ) = delete;
  
  IOmpiLevel1( IOmpiLevel1&& other ) = default;
  
  IOmpiLevel1& operator=( IOmpiLevel1&& other ) noexcept {
    IOmpiLevel1 tmp{ std::move( other ) };
    swap( other );
    return *this;
  }
  
  void swap( IOmpiLevel1& other ) noexcept {
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

void swap( IOmpiLevel1& a, IOmpiLevel1& b ) noexcept;

#endif /* IOMPILEVEL0_H_ */
