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
 *
 *   Modified: 2021-2024
 *     Author: Gregor Weiss
 *
 */
#include <mpi.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include <cmath>
#include <string_view>
#include <chrono>
#include <ctime>

#include "helper.h"
#include "HeatTransfer.h"
#include "IO.h"
#include "Settings.h"

#include "ndarray.h"

void printUsage() {
  std::cout << "Usage: heatTransfer config.xml output.bp scheme N M nx ny steps iterations"
               "iterations\n"
            << "  config: XML config file to use\n"
            << "  output: name of output data file\n"
            << "  scheme: name of IO scheme\n"
            << "  N:      number of processes in X dimension\n"
            << "  M:      number of processes in Y dimension\n"
            << "  L:      number of processes in Z dimension\n"
            << "  nx:     local array size in X dimension per processor\n"
            << "  ny:     local array size in Y dimension per processor\n"
            << "  nz:     local array size in Z dimension per processor\n"
            << "  steps:  the total number of steps to output\n"
            << "  iterations: one step consist of this many iterations\n\n";
}


void printPerf( std::string_view identifier, const double timing, const Settings& settings ) {
  auto now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::cout << identifier
            << " max. time [s] " << timing
            << " perf [GB/s] " << settings.globalGB / timing
            << " perf [GiB/s] " << settings.globalGiB / timing
            << "\n";
  std::cout << "    finished at " << std::ctime(&now_time);
}

void printTime( std::string_view identifier, const double timing ) {
  auto now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::cout << identifier
            << " time [s] " << timing
            << "\n";
  std::cout << "   finished at " << std::ctime(&now_time);
}

void checkEquality( std::vector<std::vector<double> >& input1,
                    std::vector<std::vector<double> >& input2 ) {
  bool equal = std::equal( std::begin( input1 ),
                           std::end( input1 ),
                           std::begin( input2 ) );
  if ( !equal ) {
    std::cout << "WARNING: read data is not equal to written data" << std::endl;
  }
}


void peakSignalToNoise( std::vector<double>& ref,
                        std::vector<double>& org ) {
  auto peakVal = std::numeric_limits<double>::lowest();
  double mse = 0;
  double psnr = 0;

  auto minVal = std::numeric_limits<double>::max();
  auto maxVal = std::numeric_limits<double>::lowest();

  for ( size_t i = 0; i < ref.size(); ++i ) {
    minVal = std::min(minVal, ref[i]);
    maxVal = std::max(maxVal, ref[i]);
    double sum = ref[i] - org[i];
    mse += sum * sum;
  }
  peakVal = std::max(peakVal, maxVal - minVal);
  mse /= static_cast<double>(ref.size());
  psnr = 20 * log10(peakVal/(2 * sqrt(mse)));
  
  std::cout << "MSE : " << mse << std::endl;
  std::cout << "PSNR : " << psnr << std::endl;
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
    ht.exchange(MPI_COMM_WORLD);

    for ( unsigned int t = 1; t <= settings.steps; ++t )
    {
      MPI_Barrier(MPI_COMM_WORLD);
      measTime = MPI_Wtime();

      ht.m_TIterations.clear();
      for ( unsigned int iter = 1; iter <= settings.iterations; ++iter )
      {
        ht.iterate();
        ht.exchange(MPI_COMM_WORLD);
        ht.store();
      }
      //ht.printT("humpa", MPI_COMM_WORLD);

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

      if ( settings.read ) {
        IO<IOVariant> istream( settings, MPI_COMM_WORLD );
        istream.chooseFormat( settings.format );
        std::vector<std::vector<double> > input( settings.iterations,
                                                 std::vector<double>( settings.ndx * settings.ndy * settings.ndz, -1.0 ) );

        MPI_Barrier(MPI_COMM_WORLD);
        measTime = MPI_Wtime();

        istream.read( t, input, settings, MPI_COMM_WORLD );

        MPI_Barrier(MPI_COMM_WORLD);
        measTime = MPI_Wtime() - measTime;

        checkEquality( input, ht.m_TIterations );
        peakSignalToNoise( input[0], ht.m_TIterations[0] );

/*
        std::cout << "input " <<std::endl;
        for ( auto& vec : input )
        {
          for ( auto& val : vec )
          { std::cout << val << ", ";}
        }
        std::cout << "\ninput \n" <<std::endl;

        std::cout << "ht.m_TIterations " <<std::endl;
        for ( auto& vec : ht.m_TIterations )
        {
          for ( auto& val : vec )
          { std::cout << val << ", ";}
        }
        std::cout << "\nht.m_TIteratinos \n" <<std::endl;
        */

        MPI_Reduce( &measTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if ( rank == 0 ) {
          printPerf( "Reading step " + std::to_string(t), maxTime, settings );
        }
      }

      if ( settings.remove ) {
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
    }

    if ( settings.remove ) {
      RemoveProcFolders( rank );
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
