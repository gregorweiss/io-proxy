/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * helper.cpp
 *
 * Created on: Mar 2022
 *     Author: Gregor Weiss
 *
 */

#include "helper.h"

#include <filesystem>
#include <iostream>

// Generate a folder per rank
std::string MakeProcFolders( int rank ) {
    std::filesystem::path path_name = "processor" + std::to_string(rank) + "/";
    std::filesystem::create_directory( path_name );
    return path_name.string();
}

// Remove processor folders
void RemoveProcFolders( int rank ) {
    std::filesystem::path path_name = "processor" + std::to_string(rank) + "/";
    std::filesystem::exists( path_name );
    std::filesystem::remove( path_name );
}

// Generate a file name from the outputfile string and the arguments
// default is add suffix if not already there
// if rank and step is specified, it will create a unique file name for that
// rank and step
std::string MakeFilename( const std::string &outputfile,
                          const std::string &suffix,
                          int rank /* = -1 */,
                          int step /* = -1 */ )
{
    std::string name{ outputfile };
    appendSuffix( name, suffix );

    const size_t pos{ name.size() - suffix.size() };
    insertIndexIfRequired( name, pos, rank );
    insertIndexIfRequired( name, pos, step );

    return name;
}

std::string& appendSuffix( std::string& name,
                           std::string const& suffix )
{
    if (
         ( name.size() <= suffix.size() )
         ||
         (
           ( name.size() > suffix.size() )
           &&
           ( name.find(suffix) != name.size() - suffix.size() )
         )
       )
    { name += suffix; }

    return name;
}

std::string& insertIndexIfRequired( std::string& name,
                                    size_t pos,
                                    int index )
{
    if (index >= 0)
    { name.insert( pos, "." + std::to_string(index) ); }

    return name;
}

int getRank( MPI_Comm communicator )
{
  int rank;
  MPI_Comm_rank( communicator, &rank );
  return rank;
}

int getNProcs( MPI_Comm communicator )
{
  int size;
  MPI_Comm_size( communicator, &size );
  return size;
}

