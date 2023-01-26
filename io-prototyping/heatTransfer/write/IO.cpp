/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#include "IOadios2.h"
#include "IOascii.h"
#include "IObinary.h"
#include "IOmpiLevel0.h"
#include "IOmpiLevel3.h"
#ifdef HAVE_SIONLIB
  #include "IOsion.h"
#endif
#include "IOstream.h"

using IOVariant = std::variant<IOadios2,
                               IOascii,
                               IObinary,
                               IOmpiLevel0,
                               IOmpiLevel3,
#ifdef HAVE_SIONLIB
                               IOsion,
#endif
                               IOstream>;

struct Format
{
  std::optional<IOVariant> operator()( const Settings& s, MPI_Comm comm, std::string& ioFormat ) {
    if ( ioFormat.compare( "adios2" ) == 0 )
    { return IOadios2{ s, comm }; }
    else if ( ioFormat.compare( "ascii" ) == 0 )
    { return IOascii{ s, comm }; }
    else if ( ioFormat.compare( "binary" ) == 0 )
    { return IObinary{ s, comm }; }
    else if ( ioFormat.compare( "level0" ) == 0 )
    { return IOmpiLevel0{ s, comm }; }
    else if ( ioFormat.compare( "level3" ) == 0 )
    { return IOmpiLevel3{ s, comm }; }
#ifdef HAVE_SIONLIB
    else if ( ioFormat.compare( "sion" ) == 0 )
    { return IOsion{ s, comm }; }
#endif
    else if ( ioFormat.compare( "stream" ) == 0 )
    { return IOstream{ s, comm }; }
    
    return IObinary{ s, comm };
  };
};

template<typename IOStrategy>
IO<IOStrategy>::IO( const Settings& settings, MPI_Comm communicator )
  : _settings{ settings }
  , _communicator{ communicator } {}

template<typename IOStrategy>
void IO<IOStrategy>::chooseFormat( std::string ioFormat ) {
  std::optional<IOStrategy> newFormat = Format{}( _settings, _communicator, ioFormat );
  if ( newFormat )
  { _ioFormat = std::move( *newFormat ); }
}

template<typename IOStrategy>
void IO<IOStrategy>::write( int step,
                            const HeatTransfer& ht,
                            const Settings& s,
                            MPI_Comm comm ) {
  std::visit(
    [ &step, &ht, &s, &comm ]( auto& ioFormat )
    {
      ioFormat.write( step, ht, s, comm );
    }, _ioFormat
  );
}

template<typename IOStrategy>
void IO<IOStrategy>::read( const int step,
                           std::vector<double>& buffer,
                           const Settings& s,
                           MPI_Comm comm ) {
  std::visit(
    [ &step, &buffer, &s, &comm ]( auto& ioFormat )
    {
      ioFormat.read( step, buffer, s, comm );
    }, _ioFormat
  );
}

template<typename IOStrategy>
void IO<IOStrategy>::remove( const int step ) {
  std::visit(
    [ &step ]( auto& ioFormat )
    {
      ioFormat.remove( step );
    }, _ioFormat
  );
}

