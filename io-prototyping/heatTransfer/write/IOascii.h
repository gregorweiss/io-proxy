/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOascii.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOASCII_H_
#define IOASCII_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <fstream>
#include <streambuf>

#include <vector>
#include <mpi.h>

class IOascii
{
 public:
  IOascii() = default;
  
  IOascii( const Settings& s, MPI_Comm comm );
  
  ~IOascii() = default;
  
  IOascii( IOascii const& other ) = delete;
  
  IOascii( IOascii&& other ) = default;
  
  IOascii& operator=( IOascii const& other ) = delete;
  
  IOascii& operator=( IOascii&& other ) noexcept {
    IOascii tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  };
  
  void swap( IOascii& other ) noexcept {
    using std::swap;
    swap( m_outputfilename, other.m_outputfilename );
    swap( _of, other._of );
    swap( _buf, other._buf );
  }
  
  void write( int step,
              const HeatTransfer& ht,
              const Settings& s,
              MPI_Comm comm );
  
  void read( const int step,
             std::vector<double>& buffer,
             const Settings& s,
             MPI_Comm comm );
  
  void remove( const int step );
 
 private:
  std::string m_outputfilename{};
  std::ofstream _of{};
  std::streambuf* _buf{ nullptr };
};

#endif /* IOASCII_H_ */
