## I/O Prototyping

Prototyping different I/O libraries and their possible I/O strategies based on the heat transfer example taken from the ADIOS2 repository.

The example has been made more extensible through the visitor pattern. The prototype can be used for time measurements of writing using the different approaches. It compares ADIOS2, std::fstream, SIONlib, and MPI-IO. The latter is used to implement various modes mixing collective vs independent I/O operations and contiguous vs irregular access patterns. See also "Using Advanced MPI: Modern Features of the Message-Passing Interface" from W. Gropp, T. Hoefler, R. Thakur and E. Lusk <https://dl.acm.org/doi/10.5555/2717108>

#### Acknowledgment
This application has been developed as part of the exaFOAM Project https://www.exafoam.eu, which has received funding from the European High-Performance Computing Joint Undertaking (JU) under grant agreement No 956416. The JU receives support from the European Union's Horizon 2020 research and innovation programme and France, Germany, Italy, Croatia, Spain, Greece, and Portugal.
