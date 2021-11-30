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
#include <vector>        // std::vector
#include <algorithm>    // std::shuffle

#include <iostream>

MPI_File fh;
MPI_Datatype filetype;
int count;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile + ".mpiio_write_all_indexed";

    // Set global size
    int globsize = s.gndx * s.gndy;

    // Set count of buffer, i.e. size of ht.data()
    count = s.ndx * s.ndy;

    MPI_Status status;
    int rank, size;
    MPI_Comm_rank(comm, &rank); 
    MPI_Comm_size(comm, &size); 

    //- Create an random index vectors for each process
    std::vector<int> indexes(count, 0);
    if (rank == 0)
    {
        std::vector<int> globindexes( globsize, 0 );
        // fill global unique index vector
        std::generate( globindexes.begin(), globindexes.end(), [i = 0] () mutable { return i++; } );
        // random shuffling
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle( globindexes.begin(), globindexes.end(), std::default_random_engine(seed) );
        // move and send data into local index vectors
        std::vector<int> locindexes( count, 0 );

        int n = 0; auto globIt = globindexes.begin();
        std::move( globIt, globIt + count, indexes.begin() );
        ++n; globIt += count; 
        for (
                ;
                n < size && globIt != globindexes.end();
                ++n , globIt += count
            )
        {
            std::move( globIt, globIt + count, locindexes.begin() );
            MPI_Send( locindexes.data(), count, MPI_INT, n, 0, comm ); 
        }

    } else {
         MPI_Recv( indexes.data(), count, MPI_INT, 0, 0, comm, &status);
    }
    std::sort( indexes.begin(), indexes.end() ); // monotonically increasing for file view

    // Create indexed MPI_Datatype
    MPI_Type_create_indexed_block
    (
        count,
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
        MPI_MODE_CREATE | MPI_MODE_WRONLY,
        MPI_INFO_NULL,
        &fh
    );

    // Set file view
    MPI_File_set_view
    (
        fh,
        0,
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
    MPI_File_write_all(fh, ht.data(), count, MPI_DOUBLE, MPI_STATUS_IGNORE);
}
