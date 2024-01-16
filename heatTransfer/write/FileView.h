/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileView.h
 *
 *  Created on: Mar 2023
 *      Author: Gregor Weiss
 */

#ifndef FILEVIEW_H_
#define FILEVIEW_H_

#include "HeatTransfer.h"
#include "Settings.h"

#include <functional>
#include <vector>
#include <mpi.h>

void subarray1D( const Settings& settings, MPI_Datatype& filetype );
void subarray2D( const Settings& settings, MPI_Datatype& filetype );
void subarray2D_contiguous( const Settings& settings, MPI_Datatype& filetype );
void darray1D( const Settings& settings, MPI_Datatype& filetype );
void darray2D( const Settings& settings, MPI_Datatype& filetype );
void darray2D_contiguous( const Settings& settings, MPI_Datatype& filetype );

class FileView
{
  MPI_Comm _communicator{};
  std::string _format;
  int _rank{};
  int _nprocs{};

 public:

  FileView( const Settings& settings,
            MPI_Comm communicator,
            std::function<void(const Settings&, MPI_Datatype&)> gen_filetype );

  ~FileView();

  FileView( FileView const& other );

  FileView( FileView&& other );

  FileView& operator=( FileView const& other ) noexcept {
    FileView tmp{ other };
    swap( tmp );
    return *this;
  };

  FileView& operator=( FileView&& other ) noexcept {
    FileView tmp{ std::move( other ) };
    swap( tmp );
    return *this;
  };

  void swap( FileView& other ) noexcept {
    using std::swap;
    swap( _communicator, other._communicator );
    swap( _rank, other._rank );
    swap( _nprocs, other._nprocs );
    MPI_Type_dup( other._filetype, &_filetype );
  }

  MPI_Datatype _filetype;
};

void swap( FileView& a, FileView& b ) noexcept;

#endif /* FILEVIEW_H_ */
