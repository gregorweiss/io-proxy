## Heat transfer example

This example solves a 2D Poisson equation for temperature in homogeneous media using finite differences.

### Installation

```
[CC=mpicc]
[CXX=mpicxx]
cmake .. -DADIOS2_DIR=_YOURADIOS2PATH_/lib64/cmake/adios2 [ -Dwith-sionlib=_YOURSIONLIBPATH_ ]
make
```

### Usage

Parallel execution with MPI
```
mpirun -np 24 heatTransfer heat_bp4.xml heat.bp adios2 4 3 2 1024 1024 1024 10 1
```

#### Options

```
Writer usage:  heatTransfer config.xml output scheme N M nx ny steps iterations
  config: XML config file to use
  output: name of output data file
  scheme: name of IO scheme
  N:      number of processes in X dimension
  M:      number of processes in Y dimension
  L:      number of processes in Z dimension
  nx:     local array size in X dimension per processor
  ny:     local array size in Y dimension per processor
  nz:     local array size in Z dimension per processor
  steps:  the total number of steps to output
  iterations: one step consist of this many iterations
```

`scheme` options
```
ADIOS2:
  adios2
POSIX:
  binary, binary_with_folders
MPI-IO:
  level0, level1
  level3_1Dsubarray, level3_2Dsubarray, level3_2Dsubarray_contiguous,
  level3_1Ddarray, level3_2Ddarray, level3_2Ddarray_contiguous
```

The schemes not relying on ADIOS2 do not use the XML config file. So just type "none" for the `config` argument.

