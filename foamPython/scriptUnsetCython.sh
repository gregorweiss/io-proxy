#!/bin/bash

scriptDir=$(pwd)


cythonDir=src/finiteVolume/fvSolution/
cythonFile="fvSolution"

cd $cythonDir
rm -r build "${cythonFile}."* > /dev/null 2>&1
cp -r "${cythonFile}_py.py" "${cythonFile}.py"

cd $scriptDir


cythonDir=src/finiteVolume/fvSolution/solvers/GaussSeidel
cythonFile="GaussSeidel"

cd $cythonDir
rm -r build "${cythonFile}."* > /dev/null 2>&1
cp -r "${cythonFile}_py.py" "${cythonFile}.py"

cd $scriptDir


cythonDir=src/finiteVolume/fvMatrices/fvm/fvSchemes/laplacian
cythonFile="laplacian"

cd $cythonDir
rm -r build "${cythonFile}."* > /dev/null 2>&1
cp -r "${cythonFile}_py.py" "${cythonFile}.py"

cd $scriptDir


cythonDir=src/finiteVolume/fvMatrices/fvm/fvSchemes/laplacian
cythonFile="boundaryContributions"

cd $cythonDir
rm -r build "${cythonFile}."* > /dev/null 2>&1
cp -r "${cythonFile}_py.py" "${cythonFile}.py"

cd $scriptDir


cythonDir=src/OpenFOAM/fields/volField/volVectorField
cythonFile="volVectorField"

cd $cythonDir
rm -r build "${cythonFile}."* > /dev/null 2>&1
cp -r "${cythonFile}_py.py" "${cythonFile}.py"

cd $scriptDir
