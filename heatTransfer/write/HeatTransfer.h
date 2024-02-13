/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * HeatTransfer.h
 *
 *  Created on: Feb 2017
 *      Author: Norbert Podhorszki
 */

#ifndef HEATTRANSFER_H_
#define HEATTRANSFER_H_

#include <mpi.h>

#include <memory>
#include <vector>

#include "Settings.h"

#include "ndarray.h"

class HeatTransfer
{
public:
    HeatTransfer(const Settings &settings); // Create two 2D arrays with ghost
                                            // cells to compute
    ~HeatTransfer() = default;
    void init(bool init_with_rank); // set up array values with either rank or
                                    // real demo values
    void iterate();                 // one local calculation step
    void heatEdges();               // reset the heat values at the global edge
    void exchange(MPI_Comm comm);   // send updates to neighbors

    // return a single value at index i,j. 0 <= i <= ndx+2, 0 <= j <= ndy+2
    double T(int i, int j) const { return m_TCurrent(i, j); }; // TODO: needs to be adapted for arbitrary dimensionanlity
    // return (1D) pointer to current T data without ghost cells, ndx*ndy elements
    std::vector<double> data_noghost() const;
    void store();

    void printT(std::string message,
                MPI_Comm comm) const; // debug: print local TCurrent on stdout

    std::vector<std::vector<double> > m_TIterations{};

private:
    const Settings &m_s;

    const double edgetemp = 3.0; // temperature at the edges of the global plate
    const double omega = 0.8; // weight (1-omega) in iteration

    ndarray<double> m_TCurrent;
    ndarray<double> m_TNext;
};

#endif /* HEATTRANSFER_H_ */
