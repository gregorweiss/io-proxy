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

#include <mpi.h>

MPI_File fh;
int mpiio_count;
int mpiio_rank, mpiio_size;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write_at";

    MPI_Comm_rank(comm, &mpiio_rank); 
    MPI_Comm_size(comm, &mpiio_size); 

    // Set count of buffer, i.e. size of ht.data()
    mpiio_count = s.ndx * s.ndy;

    // Open file
    MPI_File_open
    (
        comm,
        m_outputfilename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );
}

IO::~IO()
{
    MPI_File_close(&fh);
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{

    int offset = ( mpiio_rank + mpiio_size*step ) * mpiio_count * sizeof(double);
    MPI_File_write_at(fh, offset, ht.data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}
