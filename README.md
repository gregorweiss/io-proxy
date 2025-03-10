## I/O Prototyping

Prototyping different I/O libraries and their possible I/O strategies based on the heat transfer example taken from the ADIOS2 repository.

The example has been made more extensible through the visitor pattern. The prototype can be used for time measurements of writing using the different approaches. It compares ADIOS2, std::fstream, SIONlib, and MPI-IO. The latter is used to implement various modes mixing collective vs independent I/O operations and contiguous vs irregular access patterns. See also "Using Advanced MPI: Modern Features of the Message-Passing Interface" from W. Gropp, T. Hoefler, R. Thakur and E. Lusk <https://dl.acm.org/doi/10.5555/2717108>

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

#### Acknowledgment
This application has been developed as part of the exaFOAM Project https://www.exafoam.eu, which has received funding from the European High-Performance Computing Joint Undertaking (JU) under grant agreement No 956416. The JU receives support from the European Union's Horizon 2020 research and innovation programme and France, Germany, Italy, Croatia, Spain, Greece, and Portugal.
