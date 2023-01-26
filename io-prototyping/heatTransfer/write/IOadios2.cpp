/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOadios2.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#include "IOadios2.h"
#include "helper.h"

#include <filesystem>

IOadios2::IOadios2( const Settings& settings, MPI_Comm communicator )
  : _adios2Component{ std::make_unique<adios2::ADIOS>( settings.configfile, communicator ) }
  , _ioOutput{ _adios2Component->DeclareIO( "writer" ) }
  , _ioInput{ _adios2Component->DeclareIO( "reader" ) }
  , _outputVariable{ defineVariableBySettings( _outputVariable, _ioOutput, "T", settings ) }
  , _inputVariable{ defineVariable( _inputVariable, _ioInput, "T" ) }
  , _communicator{ communicator }
  , _outputfilename{ settings.outputfile }
  , _rank{ getRank( _communicator ) }
{
  _ioOutput.SetEngine( "BP4" );
  _ioInput.SetEngine( "BP4" );
}

void IOadios2::write( int step,
                      const HeatTransfer& ht,
                      const Settings& s,
                      MPI_Comm comm ) {
  std::vector<double> v = ht.data_noghost();

  _outputfilename = MakeFilename( s.outputfile, ".bp", -1, step );
  _engineWriter = _ioOutput.Open( _outputfilename, adios2::Mode::Write, _communicator );
  
  _engineWriter.BeginStep();
  _engineWriter.Put<double>( _outputVariable, v.data() );
  _engineWriter.EndStep();
  
  _engineWriter.Close();
}

void IOadios2::read( const int step,
                     std::vector<double>& buffer,
                     const Settings& s,
                     MPI_Comm comm ) {
  auto inputfilename = MakeFilename( s.outputfile, ".bp", -1, step );
  _engineReader = _ioInput.Open( inputfilename, adios2::Mode::Read, _communicator );
  _engineReader.BeginStep();

  _inputVariable = _ioInput.InquireVariable<double>( "T" );
  auto blocksInfo = _engineReader.BlocksInfo( _inputVariable, _engineReader.CurrentStep() );
  auto& info = blocksInfo[_rank];
  _inputVariable.SetBlockSelection( info.BlockID );
  _engineReader.Get( _inputVariable, buffer.data() );
  _engineReader.EndStep();

  _engineReader.Close();
}

void IOadios2::remove( const int step ) {
  if ( _rank == 0 )
    std::filesystem::remove_all( _outputfilename );
}

adios2::IO& IOadios2::declareIO( adios2::IO& ioToDeclare, std::string ioName ) {
  ioToDeclare = _adios2Component->DeclareIO( ioName );
  if ( !ioToDeclare.InConfigFile())
  {
    // if not defined by user, we can change the default settings
    // BPFile is the default engine
    ioToDeclare.SetEngine( "BP4" );
    ioToDeclare.SetParameters( {{ "num_threads", "1" }} );
    
    // ISO-POSIX file output is the default transport (called "File")
    // Passing parameters to the transport
#ifdef _WIN32
    ioToDeclare.AddTransport( "File", {{ "Library", "stdio" }} );
#else
    ioToDeclare.AddTransport( "File", {{ "Library", "posix" }} );
#endif
  }
  
  return ioToDeclare;
}

adios2::Variable<double>& IOadios2::defineVariableBySettings( adios2::Variable<double>& retVariable,
                                                              adios2::IO& ioComponent,
                                                              std::string identifier,
                                                              const Settings& settings ) {
  retVariable = ioComponent.DefineVariable<double>( identifier,
                                                    { settings.gndx, settings.gndy },   // global dimensions
                                                    { settings.offsx, settings.offsy }, // global offset
                                                    { settings.ndx, settings.ndy } );   // local size
  return retVariable;
}

adios2::Variable<double>& IOadios2::defineVariable( adios2::Variable<double>& retVariable,
                                                    adios2::IO& ioComponent,
                                                    std::string identifier ) {
  retVariable = ioComponent.DefineVariable<double>( identifier );
  return retVariable;
}

void swap( IOadios2& a, IOadios2& b ) noexcept {
  a.swap( b );
}
