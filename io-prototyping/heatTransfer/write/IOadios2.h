/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOadios2.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef IOADIOS2_H_
#define IOADIOS2_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include "helper.h"

#include <vector>
#include <mpi.h>
#include <adios2.h>
#include <memory>

class IOadios2
{
 public:
  IOadios2() = default;
  
  IOadios2( const Settings& settings, MPI_Comm communicator );
  
  ~IOadios2() = default;
  
  IOadios2( IOadios2 const& other ) = delete;
  
  IOadios2( IOadios2&& other ) = default;
  
  IOadios2& operator=( IOadios2 const& other ) = delete;
  
  IOadios2& operator=( IOadios2&& other ) noexcept {
    IOadios2 tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  }
  
  void swap( IOadios2& other ) noexcept {
    using std::swap;
    swap( _adios2Component, other._adios2Component );
    swap( _ioOutput, other._ioOutput );
    swap( _ioInput, other._ioInput );
    swap( _outputVariable, other._outputVariable );
    swap( _inputVariable, other._inputVariable );
    swap( _communicator, other._communicator );
    swap( _outputfilename, other._outputfilename );
    swap( _rank, other._rank );
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
  
  adios2::IO&
  declareIO( adios2::IO& ioToDeclare, std::string ioName );
  
  adios2::Variable<double>&
  defineVariableBySettings( adios2::Variable<double>& retVariable,
                            adios2::IO& ioComponent,
                            std::string identifier,
                            const Settings& settings );
  
  adios2::Variable<double>&
  defineVariable( adios2::Variable<double>& retVariable,
                  adios2::IO& ioComponent,
                  std::string identifier );
  
  template<typename T>
  adios2::Variable<T>&
  copyVariable( adios2::Variable<T>& retVariable,
                adios2::IO& ioComponent,
                std::string identifier,
                adios2::Variable<T> const& otherVariable ) {
    retVariable = ioComponent.DefineVariable<T>( identifier,
                                                 otherVariable.Shape(),   // global dimensions
                                                 otherVariable.Start(), // global offset
                                                 otherVariable.Count());   // local size
    return retVariable;
  }
  
  std::unique_ptr<adios2::ADIOS> _adios2Component;
  adios2::IO _ioOutput;
  adios2::IO _ioInput;
  adios2::Engine _engineWriter;
  adios2::Engine _engineReader;
  adios2::Variable<double> _outputVariable;
  adios2::Variable<double> _inputVariable;
  MPI_Comm _communicator;
  std::string _outputfilename;
  std::string _configfilename;
  int _rank;
  
};

void swap( IOadios2& a, IOadios2& b ) noexcept;

#endif /* IOADIOS2_H_ */
