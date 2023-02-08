/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * main.cpp
 *
 * Recreates heat_transfer.f90 (Fortran) ADIOS tutorial example in C++
 *
 * Created on: Feb 2017
 *     Author: Norbert Podhorszki
 * Modified on: Nov 2021 by Gregor Weiss
 *
 */
#include <mpi.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include <cmath>
#include <string_view>

#include "HeatTransfer.h"
#include "IO.h"
#include "Settings.h"

void printUsage() {
  std::cout << "Usage: heatTransfer  config   output  N  M   nx  ny   steps "
               "iterations\n"
            << "  config: XML config file to use\n"
            << "  output: name of output data file/stream\n"
            << "  N:      number of processes in X dimension\n"
            << "  M:      number of processes in Y dimension\n"
            << "  nx:     local array size in X dimension per processor\n"
            << "  ny:     local array size in Y dimension per processor\n"
            << "  steps:  the total number of steps to output\n"
            << "  iterations: one step consist of this many iterations\n\n";
}


void printPerf( std::string_view identifier, const double timing, const Settings& settings ) {
  std::cout << identifier
            << " max. time [s] " << timing
            << " perf [GB/s] " << settings.globalGB / timing
            << " perf [GiB/s] " << settings.globalGiB / timing
            << std::endl;
}

void printTime( std::string_view identifier, const double timing ) {
  std::cout << identifier
            << " time [s] " << timing
            << std::endl;
}

int main( int argc, char* argv[] ) {
  
  MPI_Init( &argc, &argv );
  
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  MPI_Comm_size(MPI_COMM_WORLD, &nproc );
  
  try
  {
    double measTime = 0.0; // individual processor timing
    double maxTime = 0.0;  // reduced maximum timing
    double totalTime = MPI_Wtime();

    Settings settings( argc, argv, rank, nproc );
    HeatTransfer ht( settings );
    IO<IOVariant> io( settings, MPI_COMM_WORLD);
    io.chooseFormat( settings.format );

    ht.init( false );
    ht.heatEdges();
    ht.exchange(MPI_COMM_WORLD);
    
    for ( unsigned int t = 1; t <= settings.steps; ++t )
    {
      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime();

      for ( unsigned int iter = 1; iter <= settings.iterations; ++iter )
      {
        ht.iterate();
        ht.exchange(MPI_COMM_WORLD);
        ht.heatEdges();
      }

      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime() - measTime;

      MPI_Reduce( &measTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if ( rank == 0 ) {
        printTime( "Calculation step " + std::to_string( t ), maxTime );
      }

      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime();
      
      io.write( t, ht, settings, MPI_COMM_WORLD);
      
      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime() - measTime;

      MPI_Reduce( &measTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      if ( rank == 0 ) {
        printPerf( "Writing step " + std::to_string(t), maxTime, settings );
      }

      IO<IOVariant> istream( settings, MPI_COMM_WORLD );
      istream.chooseFormat( settings.format );
      std::vector<double> input_buffer( settings.ndx * settings.ndy, -1.0 );

      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime();

      istream.read( t, input_buffer, settings, MPI_COMM_WORLD );

      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime() - measTime;

      bool equal = std::equal( std::begin( input_buffer ),
                               std::end( input_buffer ),
                               std::begin( ht.data_noghost() ) );
      if ( !equal ) {
        std::cout << "WARNING: read data is not equal to written data" << std::endl;
      }

      MPI_Reduce( &measTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if ( rank == 0 ) {
        printPerf( "Reading step " + std::to_string(t), maxTime, settings );
      }

      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime();

      io.remove( t );

      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime() - measTime;

      MPI_Reduce( &measTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if ( rank == 0 ) {
        printTime( "Removing step " + std::to_string( t ), maxTime );
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    totalTime = MPI_Wtime() - totalTime;

    MPI_Reduce( &totalTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      std::cout << "Total runtime = " << maxTime << "s\n";
    }
  }
  catch ( std::invalid_argument& e ) // command-line argument errors
  {
    std::cout << e.what() << std::endl;
    printUsage();
  }
  catch ( std::ios_base::failure& e ) // I/O failure (e.g. file not found)
  {
    std::cout << "I/O base exception caught\n";
    std::cout << e.what() << std::endl;
  }
  catch ( std::exception& e ) // All other exceptions
  {
    std::cout << "Exception caught\n";
    std::cout << e.what() << std::endl;
  }
  
  MPI_Finalize();
  return 0;
}
