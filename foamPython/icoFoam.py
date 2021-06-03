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
    Incompressible transient fluid flow solver based on the PISO algorithm,
    with laminar flow only
    
\*---------------------------------------------------------------------------*/
"""

import time as timeModule
import os

import functions as fn

from src.finiteVolume.fvMesh.fvMesh import *
from src.finiteVolume.fvMatrices.include import *
from src.OpenFOAM.fields.include import *
from src.finiteVolume.cfdTools.include import *

# Case setup file

caseSetup = "caseSetup.py"

exec(open(caseSetup).read())
caseSetupModifyTime = os.path.getmtime(caseSetup)


startClockTime = timeModule.perf_counter()


# Read in the mesh
mesh = fvMesh.constructFromPolyMeshFolder("constant/polyMesh")


clockTime = round(timeModule.perf_counter() - startClockTime, 2)
print("Mesh creation lasted for", clockTime, 's\n')


# Initialize fields
U = volVectorField.initialize("U", mesh, boundaryDefinition_U, internalField_U)
p = volScalarField.initialize("p", mesh, boundaryDefinition_p, internalField_p)
phi = surfaceScalarField.flux(U)


cumulativeContErr = 0.0
time = startTime
U.write(round(time, 10)) # Write 0 time
p.write(round(time, 10)) # Write 0 time
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
	
	"""---------------------------- U equation ------------------------------"""
	
	# Assemble the U matrix contributions
	ddtU = fvm.construct(U, 'ddt', ddtScheme, [deltaT])
	divUU = fvm.construct(U, 'div', divUScheme, [phi])
	divNuGradU = fvm.construct(U, 'laplacian', laplacianScheme, [nu])
	gradP = fvc.construct(U, 'grad', gradPScheme, [p])
	
	# U matrix
	UEqn = (ddtU + divUU) - divNuGradU

	# Solve UEqn
	(UEqn == (-gradP) ).solve(fvSolution_U)
	
	# Set U boundary values after solving
	U.setBoundaryValues(boundaryDefinition_U)
	
	nCorrPISO = 0
	
	# PISO loop
	while (nCorrPISO < PISO['nCorrectors']):

		"""------------------------------------------------------------------"""
		
		rAU = UEqn.rA()
		HU = UEqn.H()
		HbyA = rAU * HU

		constrainHbyA(HbyA, U)
		# print("HbyA.data.cellValues_[0][0] =", HbyA.data.cellValues_[0][0])
		# print("U.data.cellValues_[0][0] =", U.data.cellValues_[0][0])
		# print("UEqn.data.field_.data.cellValues_[0][0] =", UEqn.data.field_.data.cellValues_[0][0])
		"""-------------------------- p equation ----------------------------"""
		
		# Assemble the p matrix contributions
		divrAUGradp = fvm.construct(p, 'laplacian', laplacianScheme, [rAU])
		divrAUHU = fvc.construct(p, 'div', divHbyAScheme, [HbyA])
		
		# p matrix
		pEqn = (divrAUGradp == divrAUHU)
		
		# Set reference values for the p field	
		p.setRefValue(fvSolution_p)
		
		# Solve pEqn
		pEqn.solve(fvSolution_p)
		
		# Set p boundary values after solving
		p.setBoundaryValues(boundaryDefinition_p)

		"""------------------------------------------------------------------"""
		
		
		"""----------------------- Flux correction --------------------------"""
		
		phiHbyA = surfaceScalarField.flux(HbyA)
		
		gradPVolField = volVectorField.grad(p)
		gradPbyA = rAU * gradPVolField
		phigradPbyA = surfaceScalarField.flux(gradPbyA) # Unit = [m3/s]
		
		phi = phiHbyA - phigradPbyA # Unit = [m3/s]
		
		"""------------------------------------------------------------------"""
		
		"""--------------------- Momentum correction ------------------------"""
		
		# U = HbyA - gradPbyA
		U = HbyA - rAU * volVectorField.grad(p)
		
		U.updateBoundaryDefinition(boundaryDefinition_U)
		U.updateFieldName("U")
		U.setBoundaryValues(boundaryDefinition_U)
		
		UEqn.data.psi_ = U
		
		nCorrPISO += 1
	
	exec(continuityErrs)
	
	"""----------------------------------------------------------------------"""
	
	if ((time + 1e-11) % writeTime < 1e-10):
		
		# Write the result into a file
		U.write(round(time, 10))
		p.write(round(time, 10))
		phi.updateFieldName("phi")
		phi.write(round(time, 10))
	
	# Print out the execution time
	clockTime = round(timeModule.perf_counter() - startClockTime, 2)
	print("\nTotal execution time =", clockTime, 's\n')
	fn.printMemoryUsageInMB()

print("\nSimulation ended at time", round(time, 10))


