/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IObinary.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#include "IObinary.h"
#include "helper.h"

#include <cstdio>

#include <filesystem>

IObinary::IObinary( const Settings& s, MPI_Comm comm ) {
  if ( s.format.find("_with_folders") != std::string::npos ) {
    auto rank = getRank( MPI_COMM_WORLD );
    std::string foldername = MakeProcFolders( rank );
    m_outputfilename = foldername + s.outputfile;
  } else {
    m_outputfilename = s.outputfile;
  }
}

void IObinary::write( int step,
                      const HeatTransfer& ht,
                      const Settings& s,
                      MPI_Comm comm ) {
  auto filename = MakeFilename( m_outputfilename, ".dat", s.rank, step );
  _filestream.open( filename, std::ios_base::out );
  
  auto write_size = static_cast<std::streamsize>( s.ndx * s.ndy * sizeof( double ) );
  for ( const auto& iteration : ht.m_TIterations ) {
    _filestream.write( reinterpret_cast<const char*>(iteration.data()),
                       write_size );
  }
  
  _filestream.close();
}

void IObinary::read( const int step,
                     std::vector<std::vector<double> >& buffer,
                     const Settings& s,
                     MPI_Comm comm ) {
  auto filename = MakeFilename( m_outputfilename, ".dat", s.rank, step );
  _filestream.open( filename, std::ios_base::in );

  auto read_size = static_cast<std::streamsize>( s.ndx * s.ndy * sizeof( double ) );
  for ( auto& iteration : buffer ) {
    _filestream.read( reinterpret_cast<char*>( iteration.data() ),
                      read_size );
  }

  _filestream.close();
}

void IObinary::remove( const int step ) {
  auto rank = getRank( MPI_COMM_WORLD );
  auto filename = MakeFilename( m_outputfilename, ".dat", rank, step );
  std::filesystem::remove( filename.c_str() );
}
