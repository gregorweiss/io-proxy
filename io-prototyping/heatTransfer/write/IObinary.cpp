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

#include <iostream>
#include <cstdio>

IObinary::IObinary( const Settings& s, MPI_Comm comm ) {
  m_outputfilename = MakeFilename( s.outputfile, ".dat", s.rank );
}

void IObinary::write( int step,
                      const HeatTransfer& ht,
                      const Settings& s,
                      MPI_Comm comm ) {
  std::vector<double> v = ht.data_noghost();

  m_outputfilename = MakeFilename( s.outputfile, ".dat", s.rank, step );
  _filestream.open( m_outputfilename, std::ios_base::out );
  
  _filestream.write( reinterpret_cast<const char*>(v.data()),
                     static_cast<std::streamsize>(s.ndx * s.ndy * sizeof( double )));
  
  _filestream.close();
}

void IObinary::read( const int step,
                     std::vector<double>& buffer,
                     const Settings& s,
                     MPI_Comm comm ) {
  m_outputfilename = MakeFilename( s.outputfile, ".dat", s.rank, step );
  _filestream.open( m_outputfilename, std::ios_base::in );
  _filestream.read( reinterpret_cast<char*>(buffer.data()),
                    static_cast<std::streamsize>(s.ndx * s.ndy * sizeof( double )));
  _filestream.close();
}

void IObinary::remove( const int step ) {
  std::remove( m_outputfilename.c_str());
}
