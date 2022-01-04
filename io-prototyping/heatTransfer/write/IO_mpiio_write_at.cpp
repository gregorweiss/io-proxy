/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_mpiio.cpp
 *
 *  Created on: Nov 2021
 *      Author: Gregor Weiss
 *
 *  This refers to the example of Figure 7.3 in 'Using Advance MPI' of Gropp et al.
 */

#include "IO.h"

#include <string>
#include <cstdio>

#include <mpi.h>

MPI_File fh;
MPI_Comm m_comm;
int mpiio_count;
int mpiio_rank, mpiio_size;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write_at";

    MPI_Comm_rank(comm, &mpiio_rank); 
    MPI_Comm_size(comm, &mpiio_size); 

    // Set count of buffer, i.e. size of ht.data()
    mpiio_count = s.ndx * s.ndy;

    // Set store the communicator
    m_comm = comm;
}

IO::~IO() {}

void IO::open()
{
    // Open file
    MPI_File_open
    (
        m_comm,
        m_outputfilename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_RDWR,
        MPI_INFO_NULL,
        &fh
    );
}

void IO::close()
{
    MPI_File_close(&fh);
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{

    MPI_Offset offset = ( mpiio_rank + mpiio_size*step ) * mpiio_count * sizeof(double);
    MPI_File_write_at(fh, offset, ht.data_noghost().data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void IO::open_write_close(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{

    MPI_Offset offset = ( mpiio_rank + mpiio_size*step ) * mpiio_count * sizeof(double);
    MPI_File_write_at(fh, offset, ht.data_noghost().data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void IO::read(const int step, std::vector<double> &buffer, const Settings &s,
              MPI_Comm comm)
{

    MPI_Offset offset = ( mpiio_rank + mpiio_size*step ) * mpiio_count * sizeof(double);
    MPI_File_read_at(fh, offset, buffer.data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void IO::remove(const int step)
{
    if (mpiio_rank == 0)
        std::remove(m_outputfilename.c_str());
}
