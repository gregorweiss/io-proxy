/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOmpiLevel3.cpp
 *
 *  Created on: Mar 2022
 *      Author: Gregor Weiss
 *
 *  This refers to the example of Figure 7.2 in 'Using Advance MPI' of Gropp et al.
 */

#include "FileView.h"
#include "helper.h"

#include <cstdio>
#include <iostream> // cout

#include <mpi.h>

int greatest_common_divisor(int a, int b) {
    if (a == 0) return b;
    return greatest_common_divisor(b % a, a);
}

void subarray1D( const Settings& settings, MPI_Datatype& filetype ) {
  // Create a contiguous array
  // to avoid integer overflow in subarray
  // if the global array size is to large
  MPI_Datatype contiguous_array;
  MPI_Type_contiguous( static_cast<int>( settings.ndx * settings.ndy * settings.ndz ),
                       MPI_DOUBLE,
                       &contiguous_array );
  MPI_Type_commit( &contiguous_array );
  // Create 1D subarray based on the contiguous array
  int sizes[1] = { static_cast<int>( settings.npx * settings.npy * settings.ndz ) };
  int subsizes[1] = { 1 };
  int starts[1] = { static_cast<int>( settings.rank ) };
  MPI_Type_create_subarray( 1,
                            sizes,
                            subsizes,
                            starts,
                            MPI_ORDER_C,
                            contiguous_array,
                            &filetype );
  MPI_Type_commit( &filetype );
}

void subarray2D( const Settings& settings, MPI_Datatype& filetype ) {
  int sizes[3] = { static_cast<int>( settings.gndx ), static_cast<int>( settings.gndy ), static_cast<int>( settings.gndz ) };
  int subsizes[3] = { static_cast<int>( settings.ndx ), static_cast<int>( settings.ndy ), static_cast<int>( settings.ndz ) };
  int starts[3] = { static_cast<int>( settings.offsx ), static_cast<int>( settings.offsy ), static_cast<int>( settings.offsz ) };
  MPI_Type_create_subarray( 2,
                            sizes,
                            subsizes,
                            starts,
                            MPI_ORDER_C,
                            MPI_DOUBLE,
                            &filetype );
  MPI_Type_commit( &filetype );
}

void subarray2D_contiguous( const Settings& settings, MPI_Datatype& filetype ) {
  // Create a contiguous array
  int divisor = greatest_common_divisor( settings.ndx, settings.ndy );
  MPI_Datatype contiguous_array;
  MPI_Type_contiguous( divisor,
                       MPI_DOUBLE,
                       &contiguous_array );
  MPI_Type_commit( &contiguous_array );
  // Create 2D subarray based on the contiguous array
  int sizes[2] = { static_cast<int>( settings.gndx / divisor ),
                   static_cast<int>( settings.gndy / divisor ) };
  int subsizes[2] = { static_cast<int>( settings.ndx / divisor ),
                      static_cast<int>( settings.ndy / divisor ) };
  int starts[2] = { static_cast<int>( settings.offsx / divisor ),
                    static_cast<int>( settings.offsy / divisor ) };
  MPI_Type_create_subarray( 2,
                            sizes,
                            subsizes,
                            starts,
                            MPI_ORDER_C,
                            contiguous_array,
                            &filetype );
  MPI_Type_commit( &filetype );
  MPI_Type_free( &contiguous_array );
}


void darray1D( const Settings& settings, MPI_Datatype& filetype ) {
  // Create a contiguous array
  // to avoid integer overflow in subarray
  // if the global array size is to large
  MPI_Datatype contiguous_array;
  MPI_Type_contiguous( static_cast<int>( settings.ndx * settings.ndy * settings.ndz ),
                       MPI_DOUBLE,
                       &contiguous_array );
  MPI_Type_commit( &contiguous_array );
  // Create 1D subarray based on the contiguous array
  int gsizes[] = { static_cast<int>( settings.npx * settings.npy * settings.npz ) };
  int distribs[] = { MPI_DISTRIBUTE_BLOCK };
  int dargs[] = { MPI_DISTRIBUTE_DFLT_DARG };
  int psizes[] = { static_cast<int>( settings.npx * settings.npy * settings.npz ) };
  MPI_Type_create_darray( settings.nproc,
                          settings.rank,
                          1,
                          gsizes,
                          distribs,
                          dargs,
                          psizes,
                          MPI_ORDER_C,
                          contiguous_array,
                          &filetype );
  MPI_Type_commit( &filetype );
  MPI_Type_free( &contiguous_array );
}

void darray2D( const Settings& settings, MPI_Datatype& filetype ) {
  int gsizes[] = { static_cast<int>( settings.gndx ), static_cast<int>( settings.gndy ), static_cast<int>( settings.gndz ) };
  int distribs[] = { MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK };
  int dargs[] = { MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG };
  int psizes[] = { static_cast<int>( settings.npx ), static_cast<int>( settings.npy ), static_cast<int>( settings.npz ) };
  MPI_Type_create_darray( settings.nproc,
                          settings.rank,
                          2,
                          gsizes,
                          distribs,
                          dargs,
                          psizes,
                          MPI_ORDER_C,
                          MPI_DOUBLE,
                          &filetype );
  MPI_Type_commit( &filetype );
}

void darray2D_contiguous( const Settings& settings, MPI_Datatype& filetype ) {
  // Create a contiguous array
  int divisor = greatest_common_divisor( settings.ndx, settings.ndy );
  MPI_Datatype contiguous_array;
  MPI_Type_contiguous( divisor,
                       MPI_DOUBLE,
                       &contiguous_array );
  MPI_Type_commit( &contiguous_array );
  // Create 2D subarray based on the contiguous array
  int gsizes[] = { static_cast<int>( settings.gndx / divisor ),
                   static_cast<int>( settings.gndy / divisor ) };
  int distribs[] = { MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK };
  int dargs[] = { MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG };
  int psizes[] = { static_cast<int>( settings.npx ), static_cast<int>( settings.npy ) };
  MPI_Type_create_darray( settings.nproc,
                          settings.rank,
                          2,
                          gsizes,
                          distribs,
                          dargs,
                          psizes,
                          MPI_ORDER_C,
                          contiguous_array,
                          &filetype );
  MPI_Type_commit( &filetype );
}

FileView::FileView( const Settings& settings, MPI_Comm communicator,
                    std::function<void(const Settings&, MPI_Datatype&)> gen_filetype )
  : _communicator{ communicator }
  , _format{ settings.format }
  , _rank{ static_cast<int>( settings.rank ) }
  , _nprocs{ static_cast<int>( settings.nproc ) } {
  gen_filetype( settings, _filetype );
}

FileView::FileView( FileView const& other )
  : _communicator{ other._communicator }
  , _format{ other._format }
  , _rank{ other._rank }
  , _nprocs{ other._nprocs } {
    MPI_Type_dup( other._filetype, &_filetype );
}

FileView::FileView( FileView&& other )
  : _communicator{ std::move( other._communicator ) }
  , _format{ std::move( other._format ) }
  , _rank{ std::move( other._rank ) }
  , _nprocs{ std::move( other._nprocs ) } {
    MPI_Type_dup( other._filetype, &_filetype );
}

FileView::~FileView() noexcept {
  if ( _filetype )
    MPI_Type_free( &_filetype );
}

void swap( FileView& a, FileView& b ) noexcept {
  a.swap( b );
}

