/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Settings.h
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>

struct Settings
{
    // user arguments
    std::string configfile;
    std::string outputfile;
    std::string format;
    unsigned int npx;        // Number of processes in X dimension
    unsigned int npy;        // Number of processes in Y dimension
    unsigned int npz;        // Number of processes in Z dimension
    unsigned int ndx;        // Local array size in X dimension per process
    unsigned int ndy;        // Local array size in y dimension per process
    unsigned int ndz;        // Local array size in z dimension per process
    unsigned int steps;      // Number of output steps
    unsigned int iterations; // Number of computing iterations between steps
    bool read{ false };      // Switch to turn on re-reading
    bool remove{ false };    // Switch to turn on removal

    // calculated values from those arguments and number of processes
    unsigned int gndx; // Global array size in X dimension
    unsigned int gndy; // Global array size in Y dimension
    unsigned int gndz; // Global array size in Z dimension
    double localGB;    // Local array size in GB
    double globalGB;   // Global array size in GB
    double localGiB;    // Local array size in GiB
    double globalGiB;   // Global array size in GiB
    // X dim positions: rank 0, npx, 2npx... are in the same X position
    // Y dim positions: npx number of consecutive processes belong to one row
    // Z dim positions: npx*npy number of consecutive processes belong to one depth layer 
    // (npx
    // columns)
    unsigned int posx;  // Position of this process in X dimension
    unsigned int posy;  // Position of this process in Y dimension
    unsigned int posz;  // Position of this process in Z dimension
    unsigned int offsx; // Offset of local array in X dimension on this process
    unsigned int offsy; // Offset of local array in Y dimension on this process
    unsigned int offsz; // Offset of local array in Z dimension on this process

    int rank;           // MPI rank
    unsigned int nproc; // number of processors

    // neighbors by their MPI ranks, -1 if there is no such neighbor
    int rank_left;
    int rank_right;
    int rank_up;
    int rank_down;
    int rank_front;
    int rank_back;

    /** true: std::async Write, false (default): sync */
    bool async = false;

    Settings(int argc, char *argv[], int rank, int nproc);
};

#endif /* SETTINGS_H_ */
