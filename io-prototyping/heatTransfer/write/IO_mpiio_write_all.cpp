/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_mpiio.cpp
 *
 *  Created on: Nov 2021
 *      Author: Gregor Weiss
 */

#include "IO.h"

#include <string>

#include <mpi.h>

MPI_File fh;
MPI_Datatype filetype;
int mpiio_count;
int mpiio_rank, mpiio_size;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write_all";

    // Set global and local sizes
    int psizes[2];
    psizes[0] = s.npx;
    psizes[1] = s.npy;
    int gsizes[2];
    gsizes[0] = s.gndx;
    gsizes[1] = s.gndy;
    int globsize = s.gndx * s.gndy;

    // Set MPI distribution settings
    int distribs[2];
    distribs[0] = distribs[1] = MPI_DISTRIBUTE_BLOCK;
    int dargs[2];
    dargs[0] = dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;

    // Set count of buffer, i.e. size of ht.data()
    mpiio_count = s.ndx * s.ndy;

    // Create array MPI_Datatype
    MPI_Comm_rank(comm, &mpiio_rank);
    MPI_Comm_size(comm, &mpiio_size);
    MPI_Type_create_darray
    (
        mpiio_size,
        mpiio_rank,
        2,
        gsizes,
        distribs,
        dargs,
        psizes,
        MPI_ORDER_C,
        MPI_DOUBLE,
        &filetype
    );
    MPI_Type_commit(&filetype);
    
    // Open file
    MPI_File_open
    (
        comm,
        m_outputfilename.c_str(),
        MPI_MODE_CREATE | MPI_MODE_RDWR,
        MPI_INFO_NULL,
        &fh
    );
}

IO::~IO()
{
    MPI_File_close(&fh);
    MPI_Type_free(&filetype);
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
    // Set file view
    int offset = mpiio_size * step * mpiio_count * sizeof(double);
    MPI_File_set_view
    (
        fh,
        offset,
        MPI_DOUBLE,
        filetype,
        "native",
        MPI_INFO_NULL
    );
    MPI_File_write_all(fh, ht.data_noghost().data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void IO::read(const int step, std::vector<double> &buffer, const Settings &s,
              MPI_Comm comm)
{
    // Set file view
    int offset = mpiio_size * step * mpiio_count * sizeof(double);
    MPI_File_set_view
    (
        fh,
        offset,
        MPI_DOUBLE,
        filetype,
        "native",
        MPI_INFO_NULL
    );
    MPI_File_read_all(fh, buffer.data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}
