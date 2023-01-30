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

int main( int argc, char* argv[] ) {
  
  MPI_Init( &argc, &argv );
  
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );
  MPI_Comm_size(MPI_COMM_WORLD, &nproc );
  
  try
  {
    double timeStart = MPI_Wtime();
    Settings settings( argc, argv, rank, nproc );
    HeatTransfer ht( settings );
    IO<IOVariant> io( settings, MPI_COMM_WORLD);
    io.chooseFormat( settings.format );
    
    ht.init( false );
    ht.heatEdges();
    ht.exchange(MPI_COMM_WORLD);
    
    for ( unsigned int t = 1; t <= settings.steps; ++t )
    {
      for ( unsigned int iter = 1; iter <= settings.iterations; ++iter )
      {
        ht.iterate();
        ht.exchange(MPI_COMM_WORLD);
        ht.heatEdges();
      }
      
      double mytime = 0.0;
      auto GB = static_cast<double>(settings.gndx * settings.gndy * sizeof( double )) / 1.0e9;
      auto lGB = static_cast<double>(settings.ndx * settings.ndy * sizeof( double )) / 1.0e9;
      auto GiB = static_cast<double>(settings.gndx * settings.gndy * sizeof( double )) / std::pow( 1024.0, 3.0 );
      auto lGiB = static_cast<double>(settings.ndx * settings.ndy * sizeof( double )) / std::pow( 1024.0, 3.0 );
      
      MPI_Barrier(MPI_COMM_WORLD);
      mytime = MPI_Wtime();
      
      io.write( t, ht, settings, MPI_COMM_WORLD);
      
      MPI_Barrier(MPI_COMM_WORLD);
      mytime = MPI_Wtime() - mytime;
      
      double maxtime = 0.0;
      double mintime = 0.0;
      double avgtime = 0.0;
      MPI_Reduce( &mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce( &mytime, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce( &mytime, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      
      if ( rank == 0 )
      {
        avgtime /= nproc;
        std::cout << "Writing step " << t
                  << " max. time [s] " << maxtime
                  << " min. time [s] " << mintime
                  << " avg. time [s] " << avgtime
                  << " global size [GB] " << GB
                  << " local size [GB] " << lGB
                  << " perf [GB/s] " << GB / maxtime
                  << " global size [GiB] " << GiB
                  << " local size [GiB] " << lGiB
                  << " perf [GiB/s] " << GiB / maxtime
                  << std::endl;
      }

      // Test for content in file
      if (
           settings.format.compare( "adios2" ) == 0
           ||
           settings.format.compare( "binary" ) == 0
           ||
           settings.format.compare( "level0" ) == 0
           ||
           settings.format.compare( "level3" ) == 0
         )
      {
        IO<IOVariant> istream( settings, MPI_COMM_WORLD );
        istream.chooseFormat( settings.format );
        std::vector<double> input_buffer( settings.ndx * settings.ndy, -1.0 );

        istream.read( t, input_buffer, settings, MPI_COMM_WORLD );

        bool equal = std::equal( std::begin( input_buffer ),
                                 std::end( input_buffer ),
                                 std::begin( ht.data_noghost() ) );
        if ( !equal ) {
          std::cout << "WARNING: read data is not equal to written data" << std::endl;
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      io.remove( t );
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    //double timeEnd = MPI_Wtime();
    //if (rank == 0)
    //std::cout << "Total runtime = " << timeEnd - timeStart << "s\n";
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
