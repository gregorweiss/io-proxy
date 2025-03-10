/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * heatTransfer.cpp
 *
 * Recreates heat_transfer.f90 (Fortran) in C++
 *
 * Created on: Feb 2017
 *     Author: Norbert Podhorszki
 *
 */

#include <cstring>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdexcept>
#include <string>

//#include "ndarray.h"

#include "HeatTransfer.h"

HeatTransfer::HeatTransfer(const Settings &settings)
: m_s(settings),
  m_TCurrent{ {m_s.ndx+2, m_s.ndy+2, m_s.ndz+2 }, edgetemp },
  m_TNext{ {m_s.ndx+2, m_s.ndy+2, m_s.ndz+2 }, edgetemp }
{}

void HeatTransfer::init(bool init_with_rank)
{
    const double hx = 2.0 * 4.0 * atan(1.0) / m_s.ndx;
    const double hy = 2.0 * 4.0 * atan(1.0) / m_s.ndy;
    const double hz = 2.0 * 4.0 * atan(1.0) / m_s.ndz;

    double x, y, z;
    for (unsigned int i = 1; i <= m_s.ndx ; i++)
    {
        x = 0.0 + hx * (i - 1);
        for (unsigned int j = 1; j <= m_s.ndy ; j++)
        {
            y = 0.0 + hy * (j - 1);
            for (unsigned int k = 1 ; k <= m_s.ndz; k++)
            {
                z = 0.0 + hz * (k - 1);
                m_TCurrent(i,j,k) = 10*(cos(8 * x) + cos(6 * x) - cos(4 * x) +
                                    cos(2 * x) - cos(x) +
                                    sin(8 * y) - sin(6 * y) + sin(4 * y) -
                                    sin(2 * y) + sin(y) +
                                    cos(8 * z) + cos(6 * z) - cos(4 * z) +
                                    cos(2 * z) - cos(z));
            }
        }
    }
}

void HeatTransfer::printT(std::string message, MPI_Comm comm) const
{
    int rank, size;
    int tag = 1;
    int token;
    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank > 0)
    {
        MPI_Recv(&token, 1, MPI_INT, rank - 1, tag, comm, &status);
    }

    std::cout << "Rank " << rank << " " << message << std::endl;
    for (unsigned int i = 0; i < m_s.ndx + 2; i++)
    {
        for (unsigned int k = 0; k < m_s.ndz + 2; k++)
        {
           std::cout << "  m_TCurrent[" << i << "][][" << k << "] = ";
           for (unsigned int j = 0; j < m_s.ndy + 2; j++)
           {
               std::cout << std::setw(6) << std::setprecision(2) << m_TCurrent(i, j, k);
           }
        }
        std::cout << std::endl;
    }
    std::cout << std::flush << std::endl;

    if (rank < size - 1)
    {
        MPI_Send(&token, 1, MPI_INT, rank + 1, tag, comm);
    }
}

void HeatTransfer::iterate()
{
    for (unsigned int i = 1; i <= m_s.ndx; ++i)
    {
        for (unsigned int j = 1; j <= m_s.ndy; ++j)
        {
            for (unsigned int k = 1; k <= m_s.ndz; k++)
            {
                m_TNext(i, j, k) = omega / 6 *
                                      ( m_TCurrent(i - 1, j, k) + m_TCurrent(i + 1, j, k) +
                                        m_TCurrent(i, j - 1, k) + m_TCurrent(i, j + 1, k) +
                                        m_TCurrent(i, j, k - 1) + m_TCurrent(i, j, k + 1) ) +
                                  (1.0 - omega) * m_TCurrent(i, j, k);
            }
        }
    }
    swap(m_TCurrent, m_TNext);
}

void HeatTransfer::exchange(MPI_Comm comm)
{
    // Build a custom MPI type for the column vector to allow strided access
    MPI_Datatype tColumnVector;
    MPI_Type_vector(m_s.ndx + 2, 1, m_s.ndy + 2, MPI_REAL8, &tColumnVector);
    MPI_Type_commit(&tColumnVector);

    MPI_Datatype tColumnDepthPlain;
    MPI_Type_vector(m_s.ndx + 2, m_s.ndz + 2, (m_s.ndy + 2)*(m_s.ndz + 2), MPI_REAL8, &tColumnDepthPlain);
    MPI_Type_commit(&tColumnDepthPlain);

    // Exchange ghost cells, in the order left-right-up-down

    // send to left + receive from right
    int tag = 1;
    MPI_Status status;
    if (m_s.rank_left >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send left to rank "
        //          << m_s.rank_left << std::endl;
        MPI_Send(&m_TCurrent(0,1,0), 1, tColumnDepthPlain, m_s.rank_left, tag, comm);
    }
    if (m_s.rank_right >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from right from rank "
        //          << m_s.rank_right << std::endl;
        MPI_Recv(&m_TCurrent(0,m_s.ndy+1,0), 1, tColumnDepthPlain,
                 m_s.rank_right, tag, comm, &status);
    }

    // send to right + receive from left
    tag = 2;
    if (m_s.rank_right >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send right to rank "
        //          << m_s.rank_right << std::endl;
        MPI_Send(&m_TCurrent(0,m_s.ndy,0), 1, tColumnDepthPlain, m_s.rank_right, tag,
                 comm);
    }
    if (m_s.rank_left >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from left from rank "
        //          << m_s.rank_left << std::endl;
        MPI_Recv(&m_TCurrent(0,0,0), 1, tColumnDepthPlain, m_s.rank_left, tag, comm,
                 &status);
    }

    // Cleanup the custom column vector type
    MPI_Type_free(&tColumnVector);
    MPI_Type_free(&tColumnDepthPlain);

    MPI_Datatype tRowDepthPlain;
    MPI_Type_vector(m_s.ndz + 2, m_s.ndy + 2, m_s.ndy + 2, MPI_REAL8, &tRowDepthPlain);
    MPI_Type_commit(&tRowDepthPlain);

    // send down + receive from above
    tag = 3;
    if (m_s.rank_down >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send down to rank "
        //          << m_s.rank_down << std::endl;
        MPI_Send(&m_TCurrent(m_s.ndx,0,0), 1, tRowDepthPlain, m_s.rank_down, tag, comm);
    }
    if (m_s.rank_up >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from above from rank "
        //          << m_s.rank_up << std::endl;
        MPI_Recv(&m_TCurrent(0,0,0), 1, tRowDepthPlain, m_s.rank_up, tag, comm, &status);
    }

    // send up + receive from below
    tag = 4;
    if (m_s.rank_up >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " send up to rank " <<
        // m_s.rank_up
        //          << std::endl;
        MPI_Send(&m_TCurrent(1,0,0), 1, tRowDepthPlain, m_s.rank_up, tag, comm);
    }
    if (m_s.rank_down >= 0)
    {
        // std::cout << "Rank " << m_s.rank << " receive from below from rank "
        //          << m_s.rank_down << std::endl;
        MPI_Recv(&m_TCurrent(m_s.ndx+1,0,0), 1, tRowDepthPlain, m_s.rank_down, tag, comm, &status);
    }

    MPI_Type_free(&tRowDepthPlain);

    MPI_Datatype tColumnRowPlain;
    MPI_Type_vector((m_s.ndx + 2)*(m_s.ndy + 2), 1, (m_s.ndz + 2), MPI_REAL8, &tColumnRowPlain);
    MPI_Type_commit(&tColumnRowPlain);

    // send to front + receive from back
    tag = 5;
    if (m_s.rank_front >= 0) {
        MPI_Send(&m_TCurrent(0,0,m_s.ndz), 1, tColumnRowPlain, m_s.rank_front, tag, comm);
    }
    if (m_s.rank_back >= 0) {
        MPI_Recv(&m_TCurrent(0,0,0), 1, tColumnRowPlain, m_s.rank_back, tag, comm, &status);
    }

    // send back + receive from front
    tag = 6;
    if (m_s.rank_back >= 0) {
        MPI_Send(&m_TCurrent(0,0,1), 1, tColumnRowPlain, m_s.rank_back, tag, comm);
    }
    if (m_s.rank_front >= 0) {
        MPI_Recv(&m_TCurrent(0,0,m_s.ndz+1), 1, tColumnRowPlain, m_s.rank_front, tag, comm, &status);
    }

    MPI_Type_free(&tColumnRowPlain);

}

/* Copies the internal ndx*ndy section of the ndx+2 * ndy+2 local array
 * into a separate contiguous vector and returns it.
 * @return A vector with ndx*ndy elements
 */
std::vector<double> HeatTransfer::data_noghost() const
{
    std::vector<double> d(m_s.ndx * m_s.ndy * m_s.ndz);
    for (unsigned int i = 1; i <= m_s.ndx; ++i)
    {
        std::memcpy(&d[(i - 1) * m_s.ndy * m_s.ndz], &m_TCurrent(i,0,0),
                    m_s.ndy * m_s.ndz * sizeof(double));
    }
    return d;
}

/* Copies the iteration into a separate public vector
 * of vectors that is publicly available.
 */
void HeatTransfer::store()
{
    auto TIteration = data_noghost();
    m_TIterations.push_back( TIteration );
}
