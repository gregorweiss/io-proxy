/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_ascii.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

static std::fstream fs;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = MakeFilename(s.outputfile, ".dat", s.rank);
    fs.open( m_outputfilename, std::ios_base::in | std::ios_base::out | std::ios_base::app );
}

IO::~IO()
{
    fs.close();
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
    auto pos = static_cast<std::streamsize>(step*s.ndx*s.ndy*sizeof(double));
    fs.seekp(pos);
    fs.write(reinterpret_cast<const char*>(ht.data_noghost().data()), static_cast<std::streamsize>(s.ndx*s.ndy*sizeof(double)));
}

void IO::read(const int step, std::vector<double> &buffer, const Settings &s,
               MPI_Comm comm)
{
    auto pos = static_cast<std::streamsize>(step*s.ndx*s.ndy*sizeof(double));
    fs.seekg(pos);
    fs.read(reinterpret_cast<char*>(buffer.data()), static_cast<std::streamsize>(s.ndx*s.ndy*sizeof(double)));
}
