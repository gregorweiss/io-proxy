Example compilation

''''''''
module load gcc/10.2.0
module load openmpi/4.1.3
module load cmake

CC=mpicc CXX=mpicxx cmake .. -DADIOS2_DIR=_YOURADIOS2PATH_/lib64/cmake/adios2 [ -Dwith-sionlib=_YOURSIONLIBPATH_ ]
make

''''''''

Example run

''''''''
mpirun -np 12 ./bin/heatTransfer_write ../examples/heatTransfer/heat_file.xml  heat.bp binary 4 3 1024 1024 10 1
''''''''

examples/heatTransfer

This example solves a 2D Poisson equation for temperature in homogeneous media
using finite differences.


Example


1. Produce an output

Writer usage:  heatTransfer  config output  N  M   nx  ny   steps iterations
  config: XML config file to use
  output: name of output data file/stream
  scheme: name of IO scheme
  N:      number of processes in X dimension
  M:      number of processes in Y dimension
  nx:     local array size in X dimension per processor
  ny:     local array size in Y dimension per processor
  steps:  the total number of steps to output
  iterations: one step consist of this many iterations

The ADIOS2 executable needs an XML config file to select the Engine used for the output. The engines are: File, BP5, BP4 and HDF5, the corresponding XML config files are in the examples/heatTransfer/ directory. The "File" engine will be BP5/BP4 or HDF5 depending on the extension of the file name. 

The schemes not relying on ADIOS2 do not use XML config files, so just type "none" for the config argument.

$  mpirun -np 12 ./bin/heatTransfer_write ../examples/heatTransfer/heat_file.xml heat.bp binary 4 3 1000 1000 10 1

Scheme options

ADIOS2:
  adios2
POSIX:
  binary, binary_with_folders
MPI-IO:
  level0, level1
  level3_1Dsubarray, level3_2Dsubarray, level3_2Dsubarray_contiguous, 
  level3_1Ddarray, level3_2Ddarray, level3_2Ddarray_contiguous
