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
    Class which contains the functions needed to solve the fvMatrix, with Cython
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np
cimport numpy as np

class fvSolution():
	
	# Default solution parameters
	minIterDefault_ 	= 0
	maxIterDefault_ 	= 1000
	toleranceDefault_	= 0
	relTolDefault_		= 0
	
	"""----------------------- Other functions ------------------------------"""		
	@classmethod
	def calculateResidual(self, fvMatrix, cmpt, normFactor):
			
		mesh = fvMatrix.data.mesh_
		
		psi = fvMatrix.data.psi_.data.cellValues_[cmpt]
		
		# Total number of cells
		nCells =  mesh.geometryData.meshSize_
		# Number of internal faces
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)

		source = fvMatrix.data.source_[cmpt]
		
		lower = fvMatrix.data.lower_[cmpt]
		diag = fvMatrix.data.diag_[cmpt]
		upper = fvMatrix.data.upper_[cmpt]
		
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		# Calculate Apsi
		Apsi = self.Amul(fvMatrix, cmpt)
		
		# Calculate normalisation factor for the first iteration
		if (normFactor == None):
			normFactor = self.normFactor(Apsi, fvMatrix, cmpt)
			
		# Initial residual guess
		resArray = np.absolute(source - Apsi)
		
		# Geometric sum of residual guess magnitudes in all cells
		resSum = self.gSumMag(resArray, mesh.geometryData.V_)
		
		return resSum/normFactor, normFactor
	
	
	def checkConvergence(residual, resInit, tolerance, relTol):
		
		if (residual > tolerance and residual/resInit > relTol):
			
			return False
		
		else:
			
			return True
	
	
	def Amul(fvMatrix, int cmpt):
		
		"""------------------- Cython declarations --------------------------"""
		cdef int nCells
		cdef int nInternalFaces
		cdef int faceIndex
		cdef int ownIndex
		cdef int neiIndex
		
		cdef np.ndarray[np.float_t, ndim=1] psi
		cdef np.ndarray[np.float_t, ndim=1] lower
		cdef np.ndarray[np.float_t, ndim=1] diag
		cdef np.ndarray[np.float_t, ndim=1] upper
		cdef np.ndarray[np.int_t, ndim=1] owner
		cdef np.ndarray[np.int_t, ndim=1] neighbour
		cdef np.ndarray[np.float_t, ndim=1] resultArray
		"""------------------------------------------------------------------"""
		
		mesh = fvMatrix.data.mesh_
		
		psi = fvMatrix.data.psi_.data.cellValues_[cmpt]
		
		nCells =  mesh.geometryData.meshSize_
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		
		lower = fvMatrix.data.lower_[cmpt]
		diag = fvMatrix.data.diag_[cmpt]
		upper = fvMatrix.data.upper_[cmpt]
		
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		resultArray = np.empty(nCells, dtype = float)
		
		resultArray = diag * psi
		
		for faceIndex in range(nInternalFaces):
			
			ownIndex = owner[faceIndex]
			neiIndex = neighbour[faceIndex]
			
			resultArray[ownIndex] += upper[faceIndex] * psi[neiIndex]
			resultArray[neiIndex] += lower[faceIndex] * psi[ownIndex]
			
		return resultArray
		
	
	@classmethod
	def normFactor(self, Apsi, fvMatrix, cmpt):
		
		psi = fvMatrix.data.psi_
		
		mesh = fvMatrix.data.mesh_
		
		source = fvMatrix.data.source_
		
		
		sumA = self.sumA(fvMatrix, cmpt)
		
		sumA *= self.gAverage(psi, cmpt)
		
		return (self.gSum													   \
				   (														   \
						(np.absolute(Apsi - sumA) + np.absolute(source[cmpt] - sumA)),\
						mesh.geometryData.V_								   \
				   ) + 1.0e-15 )
	
	
	def sumA(fvMatrix, cmpt):
		
		mesh = fvMatrix.data.mesh_
		
		nCells =  mesh.geometryData.meshSize_
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)

		lower = fvMatrix.data.lower_
		diag = fvMatrix.data.diag_
		upper = fvMatrix.data.upper_
		
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		resultArray = np.empty(nCells, dtype = float)
		
		# Vectorised
		resultArray = np.copy(diag[cmpt])
		
		for faceIndex in range(nInternalFaces):
			
			ownIndex = owner[faceIndex]
			neiIndex = neighbour[faceIndex]
			
			resultArray[ownIndex] += upper[cmpt][faceIndex]
			resultArray[neiIndex] += lower[cmpt][faceIndex]
		
		return resultArray
	
	
	def gAverage(psi, cmpt):
		
		mesh = psi.data.mesh_
		
		psiValues = psi.data.cellValues_
		
		V = mesh.geometryData.V_
		
		nCells =  mesh.geometryData.meshSize_
		
		averageV = np.average(V)
		
		result = np.empty(nCells, dtype = float)
		
		# Vectorised
		result = psiValues[cmpt] * V / averageV
		
		psiAverage = np.average(result)
		
		return psiAverage
		
	
	def gSum(field, scalingField):
		
		nCells = np.size(field)
		
		scalingFieldAverage = np.average(scalingField)
		
		sum = np.sum(field * scalingField / scalingFieldAverage)
		
		return sum
	
	
	def gSumMag(field, scalingField):
		
		nCells = np.size(field)
		
		scalingFieldAverage = np.average(scalingField)
		
		sum = np.sum(np.absolute(field * scalingField) / scalingFieldAverage)
		
		return sum
		
		
		
		

"""
// ************************************************************************* //
"""
