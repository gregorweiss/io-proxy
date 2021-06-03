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
    Gauss-Seidel solver class, with Cython.
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np
cimport numpy as np

from src.finiteVolume.fvSolution.fvSolution import *

class GaussSeidel(fvSolution):
	
	@classmethod
	def solveMatrix(self, fvMatrix, fvSolutionParameters):
		
		"""------------------- Cython declarations --------------------------"""
		cdef int maxIter
		cdef int minIter
		cdef float tolerance
		cdef float relTol
		cdef int noComponents
		cdef int nCells
		cdef int cmpt
		cdef float resInit
		cdef float currentResidual
		cdef int fStart
		cdef int fEnd
		cdef int noIter
		cdef int cellIndex
		cdef float psii
		cdef int faceIndex
		
		cdef np.ndarray[np.float_t, ndim=1] psi
		cdef np.ndarray[np.float_t, ndim=1] bPrime
		cdef np.ndarray[np.float_t, ndim=1] source
		cdef np.ndarray[np.float_t, ndim=1] lower
		cdef np.ndarray[np.float_t, ndim=1] diag
		cdef np.ndarray[np.float_t, ndim=1] upper
		cdef np.ndarray[np.int_t, ndim=1] owner
		cdef np.ndarray[np.int_t, ndim=1] neighbour
		cdef np.ndarray[np.int_t, ndim=1] rowStart
		"""------------------------------------------------------------------"""
		
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

		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		rowStart = fvMatrix.data.rowStart_
		
		"""------------------------- Solving loop ---------------------------"""
		
		for cmpt in range(noComponents):
			
			psi = fvMatrix.data.psi_.data.cellValues_[cmpt]
			
			source = fvMatrix.data.source_[cmpt]
			
			lower = fvMatrix.data.lower_[cmpt]
			diag = fvMatrix.data.diag_[cmpt]
			upper = fvMatrix.data.upper_[cmpt]
			
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
						
						bPrime = np.copy(source)
						
						fEnd = rowStart[0]
						
						# Loop over the matrix rows
						for cellIndex in range(nCells):
							
							"""----------------------------------------------"""
							# # Naive implementation
							# psi[cellIndex] = source[cellIndex]
							
							# for i in range(neighbour.size):
								# if (owner[i] == cellIndex):
									# psi[cellIndex] -= upper[i] * psi[neighbour[i]]
								# elif (neighbour[i] == cellIndex):
									# psi[cellIndex] -= lower[i] * psi[owner[i]]
									
							# psi[cellIndex] /= diag[cellIndex]
							
							# psi[cellIndex] += (1.0 - alpha) * psiOld[cellIndex]
							"""----------------------------------------------"""
							# Efficient implementation
							# Start and end of the owner array for this row
							fStart = fEnd
							fEnd = rowStart[cellIndex + 1]

							# Get the accumulated neighbour side
							psii = bPrime[cellIndex]
							
							# Accumulate the owner product side
							for faceIndex in range(fStart, fEnd):
								psii -= upper[faceIndex] * psi[neighbour[faceIndex]]

							# Finish psi for this cell
							psii /= diag[cellIndex]

							# Distribute the neighbour side using psi for this cell
							for faceIndex in range(fStart, fEnd):
								bPrime[neighbour[faceIndex]] -=					   \
									lower[faceIndex] * psii
							
							psi[cellIndex] = psii
							"""----------------------------------------------"""
						
						# Increment inner iteration counter
						noIter += 1
						# Calculate the new residual
						currentResidual, normFactor = self.calculateResidual(fvMatrix, cmpt, normFactor)
			
			print("Final residual = " + str(currentResidual) + ", No Iterations " + str(noIter))



"""
// ************************************************************************* //
"""
