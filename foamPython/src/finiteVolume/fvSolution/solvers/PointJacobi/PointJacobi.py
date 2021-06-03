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
    Point Jacobi solver class.
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

from src.finiteVolume.fvSolution.fvSolution import *
	
class PointJacobi(fvSolution):
	
	@classmethod
	def solveMatrix(self, fvMatrix, fvSolutionParameters):
		
		"""------------------ Data needed for the solver --------------------"""
		maxIter = self.maxIterDefault_
		minIter = self.minIterDefault_
		tolerance = self.toleranceDefault_
		relTol = self.relTolDefault_
		
		# Read solution parameters if they are present
		try:
			maxIter = fvSolutionParameters['maxIter']
		except:
			None
		try:
			minIter = fvSolutionParameters['minIter']
		except:
			None
		try:
			tolerance = fvSolutionParameters['tolerance']
		except:
			None
		try:
			relTol = fvSolutionParameters['relTol']
		except:
			None
		
		mesh = fvMatrix.data.mesh_
		
		noComponents = fvMatrix.data.psi_.noComponents_
		
		equationName = fvMatrix.data.psi_.data.fieldName_
		
		nCells = mesh.geometryData.meshSize_
		
		psi = fvMatrix.data.psi_.data.cellValues_
		
		source = fvMatrix.data.source_
		
		lower = fvMatrix.data.lower_
		diag = fvMatrix.data.diag_
		upper = fvMatrix.data.upper_

		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		rowStart = fvMatrix.data.rowStart_
		
		"""------------------------- Solving loop ---------------------------"""
		
		for cmpt in range(noComponents):
			
			componentName = fvMatrix.data.psi_.componentsNames_[cmpt]
			
			normFactor = None
			
			resInit, normFactor = self.calculateResidual(fvMatrix, cmpt, normFactor)
			
			print("Point Jacobi: Solving for " + equationName + componentName +		   \
					", Initial residual = " + str(resInit) + ", ", end='')
			
			currentResidual = resInit
			
			if (resInit > 0.0):
			
				if ((minIter > 0) or not self.checkConvergence(currentResidual, resInit, tolerance, relTol)):
					
					noIter = 0
					
					while ((													   \
						(noIter < maxIter)										   \
						 and not self.checkConvergence(currentResidual, resInit, tolerance, relTol))
						 or noIter < minIter):
						
						psiOld = np.copy(psi[cmpt])
						
						psi[cmpt] = np.copy(source[cmpt])
							
						# Loop through the internal faces and subtract neighbour
						# contributions
						for faceIndex in range(neighbour.size):
							
							psi[cmpt][owner[faceIndex]] -= upper[cmpt][faceIndex] *\
								psiOld[neighbour[faceIndex]]
								
							psi[cmpt][neighbour[faceIndex]] -= lower[cmpt][faceIndex] *\
								psiOld[owner[faceIndex]]
						
						# psi = 1 / ap * (b - a_n * psi_n)
						psi[cmpt] /= diag[cmpt]
						
						# Increment inner iteration counter
						noIter += 1
						# Calculate the new residual
						currentResidual, normFactor = self.calculateResidual(fvMatrix, cmpt, normFactor)
			
			print("Final residual = " + str(currentResidual) + ", No Iterations " + str(noIter))


"""
// ************************************************************************* //
"""
