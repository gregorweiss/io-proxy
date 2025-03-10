/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOstream.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOSTREAM_H_
#define IOSTREAM_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <fstream>
#include <vector>
#include <mpi.h>

class IOstream
{
 public:
  IOstream() = default;
  
  IOstream( const Settings& s, MPI_Comm comm );
  
  ~IOstream();
  
  IOstream( IOstream const& other ) = delete;
  
  IOstream& operator=( IOstream const& other ) = delete;
  
  IOstream( IOstream&& other ) noexcept
    : _filename{ std::move( other._filename ) } {}
  
  IOstream& operator=( IOstream&& other ) noexcept {
    IOstream tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  }
  
  void swap( IOstream& other ) noexcept {
    using std::swap;
    swap( _filename, other._filename );
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
  std::string _filename{};
  FILE* _filestream;
};

#endif /* IOSTREAM_H_ */
