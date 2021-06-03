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
    Boundary contributions for the fvm laplacian operator, with Cython.
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np
cimport numpy as np

cdef int noComponents
cdef int nCells
cdef int nInternalFaces
cdef int cmpt
cdef int faceIndex
cdef int boundaryCellIndex
cdef int boundaryIndex
cdef int cellIndexP
cdef int cellIndexN
cdef int nFaces
cdef int startFace
cdef int gammaCmpt


class laplacianBoundaryContributions:
	
	class boundaryContributions:
		
		class scalarGamma:
		
			def fixedValue(mesh, psi, boundaryPatch, source, diag, gamma):
				
				"""----------------- Cython declarations --------------------"""
				cdef np.ndarray[np.float_t, ndim=2] psiBoundaryValues
				
				cdef np.ndarray[np.float_t, ndim=1] magSf
				cdef np.ndarray[np.float_t, ndim=2] C
				cdef np.ndarray[np.float_t, ndim=2] Cf
				cdef np.ndarray[np.int_t, ndim=1] owner
				"""----------------------------------------------------------"""
				
				boundaryName = boundaryPatch[0]
				nFaces = boundaryPatch[2]
				startFace = boundaryPatch[3]
				
				# Face area magnitudes
				magSf = mesh.geometryData.magSf_
				# Cell centres
				C = mesh.geometryData.C_
				# Face centres
				Cf = mesh.geometryData.Cf_
				# Values of psi on the boundary
				psiBoundaryValues = psi.data.boundaryValues_
				
				noComponents = psi.noComponents_
				
				owner = mesh.connectivityData.owner_
				
				# Number of internal faces
				nInternalFaces = np.size(mesh.connectivityData.neighbour_)
				
				for cmpt in range(noComponents):
				
					for faceIndex in range(startFace, startFace + nFaces):
						
						boundaryCellIndex = owner[faceIndex]
							   
						boundaryIndex = faceIndex - nInternalFaces
						
						# ap = - gamma * |Sf| / |db|
						diag[cmpt][boundaryCellIndex] += 					   	   \
							- gamma * magSf[faceIndex] / 					       \
							(np.linalg.norm(Cf[faceIndex] - C[boundaryCellIndex]))
						
						# b = - gamma * |Sf| / |db| * psi_b
						source[cmpt][boundaryCellIndex] += 				  		   \
							- gamma * magSf[faceIndex] / 						   \
							(np.linalg.norm(Cf[faceIndex] - C[boundaryCellIndex])) \
							* psiBoundaryValues[cmpt][boundaryIndex]
				
				return source, diag
				
				
			def empty(mesh, psi, boundaryPatch, source, diag, gamma):
				
				None
				
				return source, diag
				
				
			def fixedGradient(mesh, psi, boundaryPatch, source, diag, gamma):
				
				"""----------------- Cython declarations --------------------"""
				cdef np.ndarray[np.float_t, ndim=1] gradPsi_b
				
				cdef np.ndarray[np.float_t, ndim=1] magSf
				cdef np.ndarray[np.int_t, ndim=1] owner
				"""----------------------------------------------------------"""
				
				boundaryName = boundaryPatch[0]
				nFaces = boundaryPatch[2]
				startFace = boundaryPatch[3]
				
				# Face area magnitudes
				magSf = mesh.geometryData.magSf_
				
				# Number of internal faces
				nInternalFaces = np.size(mesh.connectivityData.neighbour_)
				
				# The value of the gradient
				gradPsi_b = (psi.data.boundaryDefinition_[boundaryName])[1]
				
				noComponents = psi.noComponents_
				
				owner = mesh.connectivityData.owner_
				
				for cmpt in range(noComponents):
				
					for faceIndex in range(startFace, startFace + nFaces):
						
						boundaryCellIndex = owner[faceIndex]
							   
						boundaryIndex = faceIndex - nInternalFaces

						# b = - gamma * |Sf| * gradPsi_b
						source[cmpt][boundaryCellIndex] += 				 	 	   \
							 - gamma * magSf[faceIndex] * gradPsi_b[cmpt]
							
				return source, diag
				
		class volScalarFieldGamma:

			def fixedValue(mesh, psi, boundaryPatch, source, diag, gamma):
				
				"""----------------- Cython declarations --------------------"""
				cdef np.ndarray[np.float_t, ndim=2] gammaBoundaryValues
				
				cdef np.ndarray[np.float_t, ndim=1] magSf
				cdef np.ndarray[np.float_t, ndim=2] C
				cdef np.ndarray[np.float_t, ndim=2] Cf
				cdef np.ndarray[np.int_t, ndim=1] owner
				"""----------------------------------------------------------"""
				
				gammaBoundaryValues = gamma.data.boundaryValues_
				
				boundaryName = boundaryPatch[0]
				nFaces = boundaryPatch[2]
				startFace = boundaryPatch[3]
				
				# Face area magnitudes
				magSf = mesh.geometryData.magSf_
				# Cell centres
				C = mesh.geometryData.C_
				# Face centres
				Cf = mesh.geometryData.Cf_
				# Values of psi on the boundary
				psiBoundaryValues = psi.data.boundaryValues_
				
				noComponents = psi.noComponents_
				
				owner = mesh.connectivityData.owner_
				
				# Number of internal faces
				nInternalFaces = np.size(mesh.connectivityData.neighbour_)
				
				# Gamma is a volScalarField, so only one component
				gammaCmpt = 0
				for cmpt in range(noComponents):
				
					for faceIndex in range(startFace, startFace + nFaces):
						
						boundaryCellIndex = owner[faceIndex]
							   
						boundaryIndex = faceIndex - nInternalFaces
						
						# ap = - gamma * |Sf| / |db|
						diag[cmpt][boundaryCellIndex] += 				   	   \
							- gammaBoundaryValues[gammaCmpt][boundaryIndex] *  \
							magSf[faceIndex] / 							       \
							(np.linalg.norm(Cf[faceIndex] - C[boundaryCellIndex]))
						
						# b = - gamma * |Sf| / |db| * psi_b
						source[cmpt][boundaryCellIndex] += 					   \
							- gammaBoundaryValues[gammaCmpt][boundaryIndex] *  \
							magSf[faceIndex] / 								   \
							(np.linalg.norm(Cf[faceIndex] - C[boundaryCellIndex])) \
							* psiBoundaryValues[cmpt][boundaryIndex]
				
				return source, diag
				
				
			def empty(mesh, psi, boundaryPatch, source, diag, gamma):
				
				None
				
				return source, diag
				
				
			def fixedGradient(mesh, psi, boundaryPatch, source, diag, gamma):
				
				"""----------------- Cython declarations --------------------"""
				cdef np.ndarray[np.float_t, ndim=1] gradPsi_b
				cdef np.ndarray[np.float_t, ndim=2] gammaBoundaryValues
				
#				cdef np.ndarray[np.float_t, ndim=2] source
				cdef np.ndarray[np.float_t, ndim=1] magSf
				cdef np.ndarray[np.int_t, ndim=1] owner
				"""----------------------------------------------------------"""
				
				gammaBoundaryValues = gamma.data.boundaryValues_
				
				boundaryName = boundaryPatch[0]
				nFaces = boundaryPatch[2]
				startFace = boundaryPatch[3]
				
				# Face area magnitudes
				magSf = mesh.geometryData.magSf_
				
				# Number of internal faces
				nInternalFaces = np.size(mesh.connectivityData.neighbour_)
				
				# The value of the gradient
				gradPsi_b = (psi.data.boundaryDefinition_[boundaryName])[1]
				
				noComponents = psi.noComponents_
				
				owner = mesh.connectivityData.owner_
				
				# Gamma is a volScalarField, so only one component
				gammaCmpt = 0
				for cmpt in range(noComponents):
				
					for faceIndex in range(startFace, startFace + nFaces):
						
						boundaryCellIndex = owner[faceIndex]
							   
						boundaryIndex = faceIndex - nInternalFaces

						# b = - gamma * |Sf| * gradPsi_b
						source[cmpt][boundaryCellIndex] += 			 	 	   \
							 - gammaBoundaryValues[gammaCmpt][boundaryIndex] * \
							 magSf[faceIndex] * gradPsi_b
							
				return source, diag
"""
// ************************************************************************* //
"""
