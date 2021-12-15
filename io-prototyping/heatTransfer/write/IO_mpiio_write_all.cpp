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
    const int psizes[2] = { static_cast<const int>(s.npx), static_cast<const int>(s.npy) };
    const int gsizes[2] = { static_cast<const int>(s.gndx), static_cast<const int>(s.gndy) };
    const int globsize  =   static_cast<const int>(s.gndx * s.gndy);

    // Set MPI distribution settings
    constexpr int distribs[2] = { MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK };
    constexpr int dargs[2] = { MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG } ;

    // Set count of buffer, i.e. size of ht.data()
    mpiio_count = s.ndx * s.ndy;

    // Create array MPI_Datatype
    MPI_Comm_rank(comm, &mpiio_rank);
    MPI_Comm_size(comm, &mpiio_size);
    MPI_Datatype darraytype;
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
        &darraytype
    );
    
    MPI_Aint lb     = 0;
    MPI_Aint extend = static_cast<MPI_Aint>( globsize * sizeof(double) );
    MPI_Type_create_resized(darraytype, lb, extend, &filetype);

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

    // Set file view
    MPI_Offset offset = 0;
    MPI_File_set_view
    (
        fh,
        offset,
        MPI_DOUBLE,
        filetype,
        "native",
        MPI_INFO_NULL
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
    MPI_File_write_all(fh, ht.data_noghost().data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void IO::read(const int step, std::vector<double> &buffer, const Settings &s,
              MPI_Comm comm)
{
    MPI_File_read_all(fh, buffer.data(), mpiio_count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}
