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
    npz = convertToUint("M", argv[6]);
    ndx = convertToUint("nx", argv[7]);
    ndy = convertToUint("ny", argv[8]);
    ndz = convertToUint("nz", argv[9]);
    steps = convertToUint("steps", argv[10]);
    iterations = convertToUint("iterations", argv[11]);

    if ( argc >= 13 ) {
        if ( std::string(argv[12]) == std::string("read") ) {
            read = true;
        }
        if ( std::string(argv[12]) == std::string("remove") ) {
            remove = true;
        }
    }

    if ( argc >= 14 ) {
        if ( std::string(argv[12]) == std::string("read")
             ||
             std::string(argv[13]) == std::string("read")) {
            read = true;
        }
        if ( std::string(argv[12]) == std::string("remove")
             ||
             std::string(argv[13]) == std::string("remove")) {
            remove = true;
        }
    }


    if (npx * npy * npz != this->nproc)
    {
        throw std::invalid_argument("N*M must equal the number of processes");
    }

    // calculate global array size and the local offsets in that global space
    gndx = npx * ndx;
    gndy = npy * ndy;
    gndz = npz * ndz;
    double global_bytes = static_cast<double>(gndx) * static_cast<double>(gndy) *
                          static_cast<double>(gndz) * sizeof(double);
    double local_bytes  = static_cast<double>(ndx) * static_cast<double>(ndy) *
                          static_cast<double>(ndz) * sizeof(double);
    localGB   = local_bytes  / std::pow( 10.0, 9 ) * static_cast<double>( iterations );
    globalGB  = global_bytes / std::pow( 10.0, 9 ) * static_cast<double>( iterations );
    localGiB  = local_bytes  / std::pow( 2.0, 30 ) * static_cast<double>( iterations );
    globalGiB = global_bytes / std::pow( 2.0, 30 ) * static_cast<double>( iterations );
    posx = rank % npx;
    posy = ( rank / npx ) % npy;
    posz = rank / (npx * npy);
    offsx = posx * ndx;
    offsy = posy * ndy;
    offsz = posz * ndz;

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

    if (posz == 0)
        rank_back = -1;
    else
        rank_back = rank - (npx * npy);

    if (posz == npz - 1)
        rank_front = -1;
    else
        rank_front = rank + (npx * npy);

}
