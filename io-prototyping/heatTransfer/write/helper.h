/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * helper.h
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 */

#ifndef HELPER_H_
#define HELPER_H_

#include <string>
#include <mpi.h>

std::string MakeFilename( const std::string &outputfile,
                          const std::string &suffix,
                          int rank = -1,
                          int step = -1 );

std::string& appendSuffix( std::string& name,
                           std::string const& suffix );

std::string& insertIndexIfRequired( std::string& name,
                                    size_t pos,
                                    int index );

int getRank( MPI_Comm communicator );

int getNProcs( MPI_Comm communicator );

#endif /* HELPER_H_ */ 
