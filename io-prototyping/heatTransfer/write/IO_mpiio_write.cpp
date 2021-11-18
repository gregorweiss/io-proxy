/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_mpiio.cpp
 *
 *  Created on: Nov 2021
 *      Author: Gregor Weiss
 *
 *  This refers to the example of Figure 7.2 in 'Using Advance MPI' of Gropp et al.
 */

#include "IO.h"

#include <string>

#include <mpi.h>

MPI_File fh;
int mpiio_count;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write";

    int rank, size;
    MPI_Comm_rank(comm, &rank); 
    MPI_Comm_size(comm, &size); 

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

    // Set file pointer
    MPI_File_seek
    (
        fh,
        rank*mpiio_count,
        MPI_SEEK_SET
    );
}

IO::~IO()
{
    MPI_File_close(&fh);
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
    MPI_File_write(fh, ht.data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}
