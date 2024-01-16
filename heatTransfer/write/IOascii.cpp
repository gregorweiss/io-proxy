/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOascii.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#include "IOascii.h"
#include "helper.h"

#include <iomanip>
#include <iostream>
#include <cstdio>

IOascii::IOascii( const Settings& s, MPI_Comm comm ) {
  m_outputfilename = MakeFilename( s.outputfile, ".txt", s.rank );
}

void IOascii::write( int step,
                     const HeatTransfer& ht,
                     const Settings& s,
                     MPI_Comm comm ) {
  m_outputfilename = MakeFilename( s.outputfile, ".txt", s.rank, step );
  _of.open( m_outputfilename );
  _buf = _of.rdbuf();
  
  std::ostream out( _buf );
  if ( step == 0 )
  {
    out << "rank=" << s.rank << " size=" << s.ndx << "x" << s.ndy
        << " offsets=" << s.offsx << ":" << s.offsy << " step=" << step
        << std::endl;
    out << " time   row   columns " << s.offsy << "..."
        << s.offsy + s.ndy - 1 << std::endl;
    out << "        ";
    for ( unsigned int j = 1; j <= s.ndy; ++j )
    {
      out << std::setw( 9 ) << s.offsy + j - 1;
    }
    out << "\n-------------------------------------------------------------"
           "-\n";
  } else
  {
    out << std::endl;
  }
  
  out << std::fixed;
  for ( unsigned int i = 1; i <= s.ndx; ++i )
  {
    out << std::setw( 5 ) << step << std::setw( 5 ) << s.offsx + i - 1;
    for ( unsigned int j = 1; j <= s.ndy; ++j )
    {
      out << std::setw( 9 ) << std::setprecision( 5 ) << ht.T( i, j );
    }
    out << std::endl;
  }
  
  _of.close();
}

void IOascii::read( const int step,
                    std::vector<std::vector<double> >& ht,
                    const Settings& s,
                    MPI_Comm comm ) { std::cout << "IOascii::read not implemented for ascii format." << std::endl; }

void IOascii::remove( const int step ) {
  std::remove( m_outputfilename.c_str());
}

