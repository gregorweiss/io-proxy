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
int mpiio_rank, mpiio_size;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write";

    MPI_Comm_rank(comm, &mpiio_rank); 
    MPI_Comm_size(comm, &mpiio_size); 

    // Set count of buffer, i.e. size of ht.data()
    mpiio_count = s.ndx * s.ndy;
   
    // Open file
    MPI_File_open
    (
        comm,
        m_outputfilename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_RDWR,
        MPI_INFO_NULL,
        &fh
    );

    // Set initial rank-related offset
    MPI_Offset offset = mpiio_rank * mpiio_count * sizeof(double);
    MPI_File_seek( fh, offset, MPI_SEEK_SET );
}

IO::~IO()
{
    MPI_File_close(&fh);
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
    MPI_File_write( fh, ht.data_noghost().data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE );

    // Set file pointer to next time step (Assumes sequentiell writes.)
    // Prevents MPI_Offset overflow.
    MPI_Offset offset = ( mpiio_size - 1 ) * mpiio_count * sizeof(double);
    MPI_File_seek( fh, offset, MPI_SEEK_CUR );
}

void IO::read(const int step, std::vector<double> &buffer, const Settings &s,
              MPI_Comm comm)
{
    // We avoid extending the IO interface by a file pointer reset functionality
    // by restoring the file pointer when reading the first time step.
    // Definitly violates SRP.
    if ( step == 0 ) 
    {
        // Set initial rank-related offset
        MPI_Offset offset = mpiio_rank * mpiio_count * sizeof(double);
        MPI_File_seek( fh, offset, MPI_SEEK_SET );
    }

    MPI_File_read( fh, buffer.data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE );

    // Set file pointer to next time step. Assumes sequentiell reads!
    // Prevents MPI_Offset overflow.
    MPI_Offset offset = ( mpiio_size - 1 ) * mpiio_count * sizeof(double);
    MPI_File_seek( fh, offset, MPI_SEEK_CUR );
}
