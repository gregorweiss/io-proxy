/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.h
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 *
 *    Modified: 2022-2024
 *      Author: Gregor Weiss
 */

#ifndef IO_H_
#define IO_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>

#include <variant>
#include <optional>

#include <vector>
#include <mpi.h>

template<typename IOStrategy>
class IO
{
 public:
  IO() = default;
  
  IO( const Settings& settings, MPI_Comm communicator );
  
  void chooseFormat( std::string ioFormat );
  
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
  const Settings _settings;
  MPI_Comm _communicator;
  IOStrategy _ioFormat;
};

#include "IO.cpp"

#endif /* IO_H_ */
