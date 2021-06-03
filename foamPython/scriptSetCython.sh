#!/bin/bash

scriptDir=$(pwd)




cythonDir=src/finiteVolume/fvSolution/
cythonFile="fvSolution"

cd $cythonDir
rm -r $cythonFile.py > /dev/null 2>&1
cp -r "${cythonFile}_cy.pyx" $cythonFile.pyx
python3 "setup_${cythonFile}.py" build_ext --inplace

cd $scriptDir



cythonDir=src/finiteVolume/fvSolution/solvers/GaussSeidel/
cythonFile="GaussSeidel"

cd $cythonDir
rm -r $cythonFile.py > /dev/null 2>&1
cp -r "${cythonFile}_cy.pyx" $cythonFile.pyx
python3 "setup_${cythonFile}.py" build_ext --inplace

cd $scriptDir


# These are not a big improvement:

cythonDir=src/finiteVolume/fvMatrices/fvm/fvSchemes/laplacian
cythonFile="laplacian"

cd $cythonDir
rm -r $cythonFile.py > /dev/null 2>&1
cp -r "${cythonFile}_cy.pyx" $cythonFile.pyx
python3 "setup_${cythonFile}.py" build_ext --inplace

cd $scriptDir



cythonDir=src/finiteVolume/fvMatrices/fvm/fvSchemes/laplacian
cythonFile="boundaryContributions"

cd $cythonDir
rm -r $cythonFile.py > /dev/null 2>&1
cp -r "${cythonFile}_cy.pyx" $cythonFile.pyx
python3 "setup_${cythonFile}.py" build_ext --inplace

cd $scriptDir



cythonDir=src/OpenFOAM/fields/volField/volVectorField
cythonFile="volVectorField"

cd $cythonDir
rm -r $cythonFile.py > /dev/null 2>&1
cp -r "${cythonFile}_cy.pyx" $cythonFile.pyx
python3 "setup_${cythonFile}.py" build_ext --inplace

cd $scriptDir



