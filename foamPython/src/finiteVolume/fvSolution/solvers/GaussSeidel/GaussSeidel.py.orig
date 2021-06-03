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
    Gauss-Seidel solver class, pure Python.
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

from src.finiteVolume.fvSolution.fvSolution import *

class GaussSeidel(fvSolution):
	
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
			
			print("Gauss Seidel: Solving for " + equationName + componentName +		   \
					", Initial residual = " + str(resInit) + ", ", end='')
			
			currentResidual = resInit
			
			noIter = 0
			
			if (resInit > 0.0):
			
				if ((minIter > 0) or not self.checkConvergence(currentResidual, resInit, tolerance, relTol)):
					
					while ((													   \
						(noIter < maxIter)										   \
						 and not self.checkConvergence(currentResidual, resInit, tolerance, relTol))
						 or noIter < minIter):
						
						bPrime = np.copy(source[cmpt])
						
						fEnd = rowStart[0]
						
						# Loop over the matrix rows
						for cellIndex in range(nCells):
							
							"""-----------------------------------------"""
							# # Naive implementation
							# psi[cmpt][cellIndex] = source[cmpt][cellIndex]
							
							# for i in range(neighbour.size):
								# if (owner[i] == cellIndex):
									# psi[cmpt][cellIndex] -= upper[cmpt][i] * psi[cmpt][neighbour[i]]
								# elif (neighbour[i] == cellIndex):
									# psi[cmpt][cellIndex] -= lower[cmpt][i] * psi[cmpt][owner[i]]
									
							# psi[cmpt][cellIndex] /= diag[cmpt][cellIndex]
							
							# psi[cmpt][cellIndex] += (1.0 - alpha) * psiOld[cmpt][cellIndex]
							"""-----------------------------------------"""
							 
							# Efficient implementation
							# Start and end of the owner array for this row
							fStart = fEnd
							fEnd = rowStart[cellIndex + 1]

							# Get the accumulated neighbour side
							psii = bPrime[cellIndex]
							
							# Accumulate the owner product side
							for faceIndex in range(fStart, fEnd):
								psii -= upper[cmpt][faceIndex] * psi[cmpt][neighbour[faceIndex]]

							# Finish psi for this cell
							psii /= diag[cmpt][cellIndex]

							# Distribute the neighbour side using psi for this cell
							for faceIndex in range(fStart, fEnd):
								bPrime[neighbour[faceIndex]] -=					   \
									lower[cmpt][faceIndex] * psii
							
							psi[cmpt][cellIndex] = psii
						
						# Increment inner iteration counter
						noIter += 1
						# Calculate the new residual
						currentResidual, normFactor = self.calculateResidual(fvMatrix, cmpt, normFactor)
			
			print("Final residual = " + str(currentResidual) + ", No Iterations " + str(noIter))



"""
// ************************************************************************* //
"""
