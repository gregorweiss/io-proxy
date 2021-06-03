Author: Robert Anderluh, 2021


Python3 code (used Python 3.6.9 at the time of writing) which mimics the basic classes and functionality of OpenFOAM.

The source code is present in the src/ folder.


Implemented laplacianFoam, simpleFoam and icoFoam.


scriptProfileCode runs a profiler of the code.

scriptCompareSimple copies the results to a folder with the equivalent simulation in OpenFOAM, for comparison in ParaView or similar.

scriptSetCython does a compilation of the most computationaly intensive functions (at time of writing Amul in fvSolution and Gauss Seidel solver).

scriptUnsetCython reverts the code to a standard python / no Cython implementation.


Using the code:

Examples of Allrun and Allclean scripts are in the tutorials/ directory, with examples for different solvers.

Before running Allrun, make sure you have sourced OpenFOAM (a fork which supports blockMeshDict in the system/ directory).

All of the parameters required for running the case are defined in caseSetup.py.

Meshes should be in the OpenFOAM format and in the same location as with OpenFOAM - constant/polyMesh inside the case directory.
They can be generated using blockMesh, be inserting an existing polyMesh folder from an existing OpenFOAM case, or by using other OpenFOAM mesh generation utilities.
