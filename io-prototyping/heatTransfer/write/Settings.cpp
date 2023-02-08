/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Settings.cpp
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#include "Settings.h"

#include <errno.h>

#include <cstdlib>

#include <stdexcept>

#include <cmath>

static unsigned int convertToUint(std::string varName, char *arg)
{
    char *end;
    errno = 0;
    long retval = std::strtol(arg, &end, 10);
    if (end[0] || errno == ERANGE)
    {
        throw std::invalid_argument("Invalid value given for " + varName +
                                    ": " + std::string(arg));
    }
    if (retval < 0)
    {
        throw std::invalid_argument("Negative value given for " + varName +
                                    ": " + std::string(arg));
    }
    return static_cast<unsigned int>(retval);
}

Settings::Settings(int argc, char *argv[], int rank, int nproc) : rank{rank}
{
    if (argc < 9)
    {
        throw std::invalid_argument("Not enough arguments");
    }
    this->nproc = (unsigned int)nproc;

    configfile = argv[1];
    outputfile = argv[2];
    format = argv[3];
    npx = convertToUint("N", argv[4]);
    npy = convertToUint("M", argv[5]);
    ndx = convertToUint("nx", argv[6]);
    ndy = convertToUint("ny", argv[7]);
    steps = convertToUint("steps", argv[8]);
    iterations = convertToUint("iterations", argv[9]);

    if (npx * npy != this->nproc)
    {
        throw std::invalid_argument("N*M must equal the number of processes");
    }

    // calculate global array size and the local offsets in that global space
    gndx = npx * ndx;
    gndy = npy * ndy;
    auto global_bytes = gndx * gndy * sizeof(double);
    auto local_bytes  = ndx  * ndy  * sizeof(double);
    localGB  = static_cast<double>( local_bytes  / std::pow( 10, 9 ) );
    globalGB = static_cast<double>( global_bytes / std::pow( 10, 9 ) );
    localGiB = static_cast<double>( local_bytes  / std::pow( 2, 30 ) );
    globalGiB  = static_cast<double>( global_bytes / std::pow( 2, 30 ) );
    posx = rank % npx;
    posy = rank / npx;
    offsx = posx * ndx;
    offsy = posy * ndy;

    // determine neighbors
    if (posx == 0)
        rank_up = -1;
    else
        rank_up = rank - 1;

    if (posx == npx - 1)
        rank_down = -1;
    else
        rank_down = rank + 1;

    if (posy == 0)
        rank_left = -1;
    else
        rank_left = rank - npx;

    if (posy == npy - 1)
        rank_right = -1;
    else
        rank_right = rank + npx;
}
