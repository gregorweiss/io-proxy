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
    Class which contains the laplacian contributions for fvm, with Cython.
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np
cimport numpy as np

from src.OpenFOAM.fields.volField.volScalarField.volScalarField import *
from src.finiteVolume.fvMatrices.fvm.fvSchemes.laplacian.boundaryContributions import *

cdef int noComponents
cdef int nCells
cdef int nInternalFaces
cdef int cmpt
cdef int faceIndex
cdef int cellIndexP
cdef int cellIndexN
cdef float absPf
cdef float absNf
cdef float f
cdef float gammaInterpF
cdef int gammaCmpt

class laplacianClass:
	
	class laplacian(laplacianBoundaryContributions):
		
		@classmethod
		def linearOrthogonal(self, mesh, psi, fvVariables):
			
			gamma = fvVariables[0]
			
			# Check for type of diffusivity and call the appropriate function
			if (type(gamma) == int or type(gamma) == float):
				
				source, lower, diag, upper = self.linearOrthogonalScalar(mesh, psi, fvVariables)
				
			elif (type(gamma) == volScalarField):
				
				source, lower, diag, upper = self.linearOrthogonalVolScalarField(mesh, psi, fvVariables)
				
			else:
				
				raise RuntimeError(											   \
				"Laplacian with a diffusivity of type " + str(type(gamma)) +   \
				" has not been implemented!")
			
			return source, lower, diag, upper
			
			
		@classmethod
		def linearOrthogonalScalar(self, mesh, psi, fvVariables):
			
			"""------------------- Cython declarations ----------------------"""
			cdef float gamma
			
			cdef np.ndarray[np.float_t, ndim=1] magSf
			cdef np.ndarray[np.float_t, ndim=2] C
			cdef np.ndarray[np.float_t, ndim=2] source
			cdef np.ndarray[np.float_t, ndim=2] lower
			cdef np.ndarray[np.float_t, ndim=2] diag
			cdef np.ndarray[np.float_t, ndim=2] upper
			cdef np.ndarray[np.int_t, ndim=1] owner
			cdef np.ndarray[np.int_t, ndim=1] neighbour
			"""--------------------------------------------------------------"""
			
			gamma = fvVariables[0]
			magSf = mesh.geometryData.magSf_
			C = mesh.geometryData.C_
			
			noComponents = psi.noComponents_
			
			# Total number of cells
			nCells = mesh.read("geometryData.meshSize")
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			source = np.zeros((noComponents, nCells), dtype = float)
			
			lower = np.zeros((noComponents, nInternalFaces), dtype = float)
			diag = np.zeros((noComponents, nCells), dtype = float)
			upper = np.zeros((noComponents, nInternalFaces), dtype = float)
			
			owner = mesh.connectivityData.owner_
			neighbour = mesh.connectivityData.neighbour_
					 
			"""A MATRIX DEFINITION"""
			# Add internal faces contributions
			for cmpt in range(noComponents):
				
				for faceIndex in range(nInternalFaces):
				
					# aP = - sum_N gamma * |Sf| / |df|
					# aN = gamma * |Sf| / |df|
					
					cellIndexP = owner[faceIndex]
					cellIndexN = neighbour[faceIndex]
					
					lower[cmpt][faceIndex] = gamma * magSf[faceIndex] 		   \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
							   
					diag[cmpt][cellIndexP] -= gamma * magSf[faceIndex]		   \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
										
					diag[cmpt][cellIndexN] -= gamma * magSf[faceIndex]		   \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
										
					upper[cmpt][faceIndex] = gamma * magSf[faceIndex]		   \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
			
			
			# Add boundary contributions. The loop goes through the boundary
			# patches and, depending on the boundaryType, defines the source and
			# diagonal contributions by calling the appropriate functions
			for boundaryPatch in mesh.connectivityData.boundary_:
				boundaryName = boundaryPatch[0]
				boundaryType = (psi.data.boundaryDefinition_[boundaryName])[0]
				
				boundaryContributionsFunctionString =						   \
					"self.boundaryContributions.scalarGamma."				   \
					+ boundaryType											   \
					+ "(mesh, psi, boundaryPatch, source, diag, gamma)"
				
				source, diag = eval(boundaryContributionsFunctionString)
					
			return source, lower, diag, upper
			
		@classmethod
		def linearOrthogonalVolScalarField(self, mesh, psi, fvVariables):
			
			"""------------------- Cython declarations ----------------------"""
			cdef np.ndarray[np.float_t, ndim=2] gammaValues
			cdef np.ndarray[np.float_t, ndim=1] magSf
			cdef np.ndarray[np.float_t, ndim=2] C
			cdef np.ndarray[np.float_t, ndim=2] Cf
			cdef np.ndarray[np.float_t, ndim=2] source
			cdef np.ndarray[np.float_t, ndim=2] lower
			cdef np.ndarray[np.float_t, ndim=2] diag
			cdef np.ndarray[np.float_t, ndim=2] upper
			cdef np.ndarray[np.int_t, ndim=1] owner
			cdef np.ndarray[np.int_t, ndim=1] neighbour
			"""--------------------------------------------------------------"""
			
			gamma = fvVariables[0]
			gammaValues = gamma.data.cellValues_
			
			magSf = mesh.geometryData.magSf_
			C = mesh.geometryData.C_
			Cf = mesh.geometryData.Cf_
			
			noComponents = psi.noComponents_
			
			# Total number of cells
			nCells = mesh.read("geometryData.meshSize")
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			source = np.zeros((noComponents, nCells), dtype = float)
			
			lower = np.zeros((noComponents, nInternalFaces), dtype = float)
			diag = np.zeros((noComponents, nCells), dtype = float)
			upper = np.zeros((noComponents, nInternalFaces), dtype = float)
			
			owner = mesh.connectivityData.owner_
			neighbour = mesh.connectivityData.neighbour_
			
			# Gamma is a volScalarField, so only one component
			gammaCmpt = 0
			
			"""A MATRIX DEFINITION"""
			# Add internal faces contributions
			for cmpt in range(noComponents):
				
				for faceIndex in range(nInternalFaces):
					
					cellIndexP = owner[faceIndex]
					cellIndexN = neighbour[faceIndex]
					
					absPf = np.linalg.norm(Cf[faceIndex] - C[cellIndexP])
					absNf = np.linalg.norm(Cf[faceIndex] - C[cellIndexN])
					
					f = absNf / (absPf + absNf)

					gammaInterpF = 											   \
						f * gammaValues[gammaCmpt][cellIndexP] + 			   \
						(1 - f) * gammaValues[gammaCmpt][cellIndexN]
				
					# aP = - sum_N gamma * |Sf| / |df|
					# aN = gamma * |Sf| / |df|
					
					lower[cmpt][faceIndex] = gammaInterpF * magSf[faceIndex]   \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
							   
					diag[cmpt][cellIndexP] -= gammaInterpF * magSf[faceIndex]  \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
										
					diag[cmpt][cellIndexN] -= gammaInterpF * magSf[faceIndex]  \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
										
					upper[cmpt][faceIndex] = gammaInterpF * magSf[faceIndex]   \
								/ (np.linalg.norm(C[cellIndexN] - C[cellIndexP]))
			
			
			# Add boundary contributions. The loop goes through the boundary
			# patches and, depending on the boundaryType, defines the source and
			# diagonal contributions by calling the appropriate functions
			for boundaryPatch in mesh.connectivityData.boundary_:
				boundaryName = boundaryPatch[0]
				boundaryType = (psi.data.boundaryDefinition_[boundaryName])[0]
				
				boundaryContributionsFunctionString =						   \
					"self.boundaryContributions.volScalarFieldGamma."		   \
					+ boundaryType											   \
					+ "(mesh, psi, boundaryPatch, source, diag, gamma)"
				
				source, diag = eval(boundaryContributionsFunctionString)
					
			return source, lower, diag, upper


"""
// ************************************************************************* //
"""
