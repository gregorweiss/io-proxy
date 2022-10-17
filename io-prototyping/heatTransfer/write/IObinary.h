/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IObinary.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOBINARY_H_
#define IOBINARY_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <fstream>
#include <vector>
#include <mpi.h>

class IObinary
{
 public:
  IObinary() = default;
  
  IObinary( const Settings& s, MPI_Comm comm );
  
  ~IObinary() = default;
  
  IObinary( IObinary const& other ) = delete;
  
  IObinary& operator=( IObinary const& other ) = delete;
  
  IObinary( IObinary&& other ) = default;
  
  IObinary& operator=( IObinary&& other ) noexcept {
    IObinary tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  }
  
  void swap( IObinary& other ) noexcept {
    using std::swap;
    swap( m_outputfilename, other.m_outputfilename );
    swap( _filestream, other._filestream );
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
  std::fstream _filestream{};
  std::string m_outputfilename{};
};

#endif /* IOBINARY_H_ */
