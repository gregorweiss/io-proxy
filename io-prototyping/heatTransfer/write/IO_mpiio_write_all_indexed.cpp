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

#include <chrono>       // std::chrono::system_clock
#include <random>       // std::default_random_engine
#include <algorithm>    // std::shuffle

#include <iostream>

MPI_File fh;
MPI_Datatype filetype;
int mpiio_count;
int mpiio_rank, mpiio_size;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write_all_indexed";

    // Set global size
    const int globsize = s.gndx * s.gndy;

    // Set count of buffer, i.e. size of ht.data()
    mpiio_count = s.ndx * s.ndy;

    MPI_Status status;
    MPI_Comm_rank(comm, &mpiio_rank);
    MPI_Comm_size(comm, &mpiio_size);

    //- Create an random index vectors for each process
    std::vector<int> indexes(mpiio_count, 0);
    if (mpiio_rank == 0)
    {
        std::vector<int> globindexes( globsize, 0 );
        // fill global unique index vector
        std::generate( globindexes.begin(), globindexes.end(), [i = 0] () mutable { return i++; } );
        // random shuffling
        //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        //std::shuffle( globindexes.begin(), globindexes.end(), std::default_random_engine(seed) );
        std::random_shuffle( globindexes.begin(), globindexes.end() );
        // move and send data into local index vectors
        std::vector<int> locindexes( mpiio_count, 0 );

        int n = 0; auto globIt = globindexes.begin();
        std::move( globIt, globIt + mpiio_count, indexes.begin() );
        ++n; globIt += mpiio_count;
        for (
                ;
                n < mpiio_size && globIt != globindexes.end();
                ++n , globIt += mpiio_count
            )
        {
            std::move( globIt, globIt + mpiio_count, locindexes.begin() );
            MPI_Send( locindexes.data(), mpiio_count, MPI_INT, n, 0, comm );
        }

    } else {
         MPI_Recv( indexes.data(), mpiio_count, MPI_INT, 0, 0, comm, &status);
    }
    std::sort( indexes.begin(), indexes.end() ); // monotonically increasing for file view

    // Create indexed MPI_Datatype
    MPI_Type_create_indexed_block
    (
        mpiio_count,
        1,
        indexes.data(),
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
    MPI_Offset offset = mpiio_size * step * mpiio_count * sizeof(double);
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
    MPI_Offset offset = mpiio_size * step * mpiio_count * sizeof(double);
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
