"""
/*---------------------------------------------------------------------------*\
 ####                 ####     #   #             |
 #                    #  #     ##  #             | FoamPython
 ##  #### ####  ##### #### # # #   #### #### ### | v1.0
 #   #  # #  #  # # # #     #  # # #  # #  # # # |
 #   #### ##### # # # #    #    #  #  # #### # # |
-------------------------------------------------------------------------------

Author
	Robert Anderluh, 2021
	
Description
    Transient Laplace's equation solver
    
\*---------------------------------------------------------------------------*/
"""

import time as timeModule
import os

import functions as fn

from src.finiteVolume.fvMesh.fvMesh import *
from src.finiteVolume.fvMatrices.include import *
from src.OpenFOAM.fields.include import *

# Case setup file

caseSetup = "caseSetup.py"

exec(open(caseSetup).read())
caseSetupModifyTime = os.path.getmtime(caseSetup)

startClockTime = timeModule.perf_counter()


# Read in the mesh
mesh = fvMesh.constructFromPolyMeshFolder("constant/polyMesh")


clockTime = round(timeModule.perf_counter() - startClockTime, 2)
print("Mesh creation lasted for", clockTime, 's\n')


# Initialize field T
T = volScalarField.initialize("T", mesh, boundaryDefinition_T, internalField_T)


time = startTime
T.write(round(time, 10)) # Write 0 time
while (endTime - time > 1e-10):
	
	# Reread caseSetup if it is modified
	if ((os.path.getmtime(caseSetup) != caseSetupModifyTime) and		   \
		runTimeModifiable):
		
		print("Rereading case setup file.")
		exec(open(caseSetup).read())
		caseSetupModifyTime = os.path.getmtime(caseSetup)
		
	# Increment time
	if (time + deltaT > endTime):
		deltaT = endTime - time
	time += deltaT
	
	print("\nTime =", round(time, 10))
	
	# Assemble new ddt matrix contributions
	ddtMatrix = fvm.construct(T, 'ddt', ddtScheme, [deltaT])
	# Assemble the Laplacian matrix contributions
	laplacianMatrix = fvm.construct(T, 'laplacian', laplacianScheme, [DT])
	
	# The actual matrix to solve for
	matrix = ddtMatrix - laplacianMatrix
	
	matrix.solve(fvSolution_T)
	
	gradT = volVectorField.grad(T)
	
	if ((time + 1e-11) % writeTime < 1e-10):
		
		# Write the result into a file
		T.write(round(time, 10))
		gradT.write(round(time, 10))
	
	# Print out the execution time
	clockTime = round(timeModule.perf_counter() - startClockTime, 2)
	print("\nTotal execution time =", clockTime, 's\n')

print("\nSimulation ended at time", round(time, 10))


