/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO_ADIOS2.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "IO.h"

#include <string>

#include <iostream> // std::cout

#include <adios2.h>

adios2::ADIOS ad;
adios2::Engine bpWriter;
adios2::Engine bpReader;
adios2::Variable<double> varTout;
adios2::Variable<double> varTinp;
adios2::Variable<unsigned int> varGndx;
int mpiio_rank;

IO::IO(const Settings &s, MPI_Comm comm)
{
    m_outputfilename = s.outputfile;

    ad = adios2::ADIOS(s.configfile, comm);

    MPI_Comm_rank(comm, &mpiio_rank);

    adios2::IO bpout = ad.DeclareIO("writer");
    if (!bpout.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        bpout.SetEngine("BPFile");
        bpout.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
#ifdef _WIN32
        bpout.AddTransport("File", {{"Library", "stdio"}});
#else
        bpout.AddTransport("File", {{"Library", "posix"}});
#endif
    }
    adios2::IO bpinp = ad.DeclareIO("reader");
    if (!bpinp.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        bpinp.SetEngine("BPFile");
        bpinp.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
#ifdef _WIN32
        bpinp.AddTransport("File", {{"Library", "stdio"}});
#else
        bpinp.AddTransport("File", {{"Library", "posix"}});
#endif
    }

    // define T as 2D global array
    varTout = bpout.DefineVariable<double>
    (
        "T",
        // Global dimensions
        {s.gndx, s.gndy},
        // starting offset of the local array in the global space
        {s.offsx, s.offsy},
        // local size, could be defined later using SetSelection()
        {s.ndx, s.ndy}
    );
    varTinp = bpinp.DefineVariable<double>( "T" );

    //if (bpio.EngineType() == "BP3")
    //{
        //varT.SetMemorySelection({{1, 1}, {s.ndx + 2, s.ndy + 2}});
    //}

    bpWriter = bpout.Open(m_outputfilename, adios2::Mode::Write, comm);
    bpReader = bpinp.Open(m_outputfilename, adios2::Mode::Read, comm);

    // Promise that we are not going to change the variable sizes nor add new
    // variables
    bpWriter.LockWriterDefinitions();
    bpReader.LockReaderSelections();
}

IO::~IO()
{
    bpWriter.Close();
    bpReader.Close();
}

void IO::write(int step, const HeatTransfer &ht, const Settings &s,
               MPI_Comm comm)
{
    // using PutDeferred() you promise the pointer to the data will be intact
    // until the end of the output step.
    // We need to have the vector object here not to destruct here until the end
    // of function.
    // added support for MemorySelection
    if (bpWriter.Type() == "BP3")
    {
        bpWriter.BeginStep();
        bpWriter.Put<double>(varTout, ht.data_noghost().data()); //ht.data());
        bpWriter.EndStep();
    }
    else
    {
        bpWriter.BeginStep();
        std::vector<double> v = ht.data_noghost();
        bpWriter.Put<double>(varTout, v.data());
        bpWriter.EndStep();
    }
}

void IO::read(const int step, std::vector<double> &buffer, const Settings &s,
               MPI_Comm comm)
{
    bpReader.BeginStep( adios2::StepMode::Read );
    auto blocksInfo = bpReader.BlocksInfo( varTinp, bpReader.CurrentStep() );
    auto& info = blocksInfo[ mpiio_rank ];
    varTinp.SetBlockSelection( info.BlockID );
    bpReader.Get( varTinp, buffer.data() );
    bpReader.EndStep();
}
